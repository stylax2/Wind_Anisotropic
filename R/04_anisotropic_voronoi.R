##############################################################################
# 04_anisotropic_voronoi.R
#
# Anisotropic cost-distance Voronoi allocation:
#   - Build conductance transition matrix from friction surface
#   - Apply directional penalty based on prevailing wind (cosine weight)
#   - Competitive allocation: assign each cell to nearest station
#   - Repeat independently for 4 seasons
#
# Input:  01_intermediate/friction_*.tif (4 seasonal friction surfaces)
#         01_intermediate/station_seasonal_wind.gpkg
# Output: 00_data/vector/voronoi_aniso_spring.gpkg ... _winter.gpkg
#
# Note: This script requires significant RAM (>16 GB) and CPU time.
#       Uses doParallel with chunked processing.
#       If outputs already exist in 00_data/vector/, skip this script.
#
# References:
#   gdistance (van Etten, 2017) - cost-distance computation
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  04. Anisotropic Cost-Distance Voronoi Allocation\n")
cat("================================================================\n\n")

# --- Check if outputs already exist ---
existing <- sapply(SEASONS, function(s) {
  file.exists(vector_path(sprintf("voronoi_aniso_%s.gpkg", tolower(s))))
})

if (all(existing)) {
  cat("[SKIP] All 4 seasonal Voronoi files already exist.\n")
  cat("  Delete them from 00_data/vector/ to force re-computation.\n")
  cat("  Proceeding to next script.\n")
} else {
  
  library(raster)      # Required for gdistance compatibility
  library(gdistance)
  library(Matrix)
  library(doParallel)
  library(foreach)
  
  # =============================================================================
  # 1. Load station data
  # =============================================================================
  cat("[1/3] Loading station data...\n")
  
  stations <- st_read(inter_path("station_seasonal_wind.gpkg"), quiet = TRUE)
  coords   <- st_coordinates(stations)
  n_st     <- nrow(coords)
  
  cat(sprintf("  Stations: %d\n", n_st))
  
  # Parallel setup
  N_CORES    <- min(12, parallel::detectCores() - 2)
  CHUNK_SIZE <- 25
  chunks     <- split(seq_len(n_st), ceiling(seq_len(n_st) / CHUNK_SIZE))
  
  cl_par <- makeCluster(N_CORES)
  registerDoParallel(cl_par)
  cat(sprintf("  Parallel cores: %d | Chunks: %d\n", N_CORES, length(chunks)))
  
  # =============================================================================
  # 2. Seasonal loop
  # =============================================================================
  for (current_season in SEASONS) {
    
    out_file <- vector_path(sprintf("voronoi_aniso_%s.gpkg", tolower(current_season)))
    
    if (file.exists(out_file)) {
      cat(sprintf("\n[SKIP] %s already exists.\n", basename(out_file)))
      next
    }
    
    cat(sprintf("\n[2/3] Processing %s...\n", current_season))
    
    # Load friction surface
    cost_terra <- rast(inter_path(sprintf("friction_%s.tif", tolower(current_season))))
    cost_r     <- raster(cost_terra)
    n_total    <- ncell(cost_r)
    wd_col     <- paste0("wd_", current_season)
    
    # Initialize competition rasters
    min_cost <- setValues(cost_r, Inf)
    min_idx  <- setValues(cost_r, 0L)
    
    # Base transition matrix (conductance = 1/cost)
    cond_layer <- 1 / cost_r
    tr_base    <- transition(cond_layer, transitionFunction = mean, directions = 8)
    tr_base    <- geoCorrection(tr_base, type = "c")
    
    # Pre-compute adjacency and flow directions
    adj <- adjacent(cost_r, cells = 1:n_total, pairs = TRUE, directions = 8)
    i_idx <- adj[, 1]
    j_idx <- adj[, 2]
    
    coords_from   <- xyFromCell(cost_r, i_idx)
    coords_to     <- xyFromCell(cost_r, j_idx)
    flow_rad_vec  <- atan2(coords_to[,1] - coords_from[,1],
                           coords_to[,2] - coords_from[,2])
    flow_rad_vec[flow_rad_vec < 0] <- flow_rad_vec[flow_rad_vec < 0] + 2 * pi
    
    base_values_vec <- transitionMatrix(tr_base)[adj]
    
    # Chunked parallel processing
    for (ch in seq_along(chunks)) {
      idx_chunk <- chunks[[ch]]
      
      cost_list <- foreach(j = seq_along(idx_chunk),
                           .packages = c("gdistance", "raster", "Matrix")) %dopar% {
        i <- idx_chunk[j]
        dom_wd <- stations[[wd_col]][i]
        upwind_rad <- ((dom_wd + 180) %% 360) * (pi / 180)
        
        angle_diff <- abs(flow_rad_vec - upwind_rad)
        mask_pi <- angle_diff > pi
        angle_diff[mask_pi] <- 2 * pi - angle_diff[mask_pi]
        
        penalty_weight <- (1 + cos(angle_diff)) / 2
        penalty_weight[penalty_weight < 0.05] <- 0.05
        
        new_sparse <- Matrix::sparseMatrix(
          i = i_idx, j = j_idx,
          x = base_values_vec * penalty_weight,
          dims = c(n_total, n_total)
        )
        
        tr_aniso <- tr_base
        transitionMatrix(tr_aniso) <- new_sparse
        
        cost_i <- accCost(tr_aniso, coords[i, , drop = FALSE])
        return(cost_i[])
      }
      
      # Competitive allocation
      for (j in seq_along(idx_chunk)) {
        i <- idx_chunk[j]
        cost_vals <- cost_list[[j]]
        mask <- cost_vals < min_cost[]
        mask[is.na(mask)] <- FALSE
        min_cost[mask] <- cost_vals[mask]
        min_idx[mask]  <- i
      }
      
      rm(cost_list); gc()
      cat(sprintf("   Chunk %d/%d done\n", ch, length(chunks)))
    }
    
    # Convert to polygons
    min_idx_terra <- rast(min_idx)
    crs(min_idx_terra) <- crs(cost_terra)
    
    voronoi_poly <- as.polygons(min_idx_terra, dissolve = TRUE) |>
      st_as_sf() |>
      dplyr::rename(station_idx = 1) |>
      dplyr::mutate(station_idx = as.integer(station_idx)) |>
      dplyr::filter(!is.na(station_idx), station_idx >= 1, station_idx <= n_st) |>
      dplyr::mutate(
        site   = stations$site[station_idx],
        season = current_season
      )
    
    st_write(voronoi_poly, out_file, delete_layer = TRUE, quiet = TRUE)
    cat(sprintf("  → %s (%d polygons)\n", basename(out_file), nrow(voronoi_poly)))
  }
  
  stopCluster(cl_par)
}

cat("\n[Complete] Anisotropic Voronoi allocation done.\n")
