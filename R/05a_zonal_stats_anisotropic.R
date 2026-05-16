##############################################################################
# 05a_zonal_stats_anisotropic.R
#
# Extract zonal statistics from anisotropic Voronoi polygons:
#   - Continuous variables: mean, stdev, max, min (8 layers × 4 stats = 32)
#   - Categorical variables: area fraction (FRTP 15cls + TPI 6cls × 3 = 33)
#   - Derived: DEM_Range = max_DEM − min_DEM
#   - Season as dummy variable
#   - Filter to 206 common stations
#
# Input:  00_data/vector/voronoi_aniso_*.gpkg (4 seasons)
#         00_data/raster/*.tif (11 layers)
# Output: 01_intermediate/spatial_features_aniso.rds
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  05a. Zonal Statistics — Anisotropic Voronoi\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load rasters and build stacks
# =============================================================================
cat("[1/4] Loading rasters...\n")

dem   <- rast(raster_path("dem.tif"));        names(dem)   <- "DEM"
slope <- rast(raster_path("slope.tif"));      names(slope) <- "Slope"

aspect     <- rast(raster_path("aspect.tif"))
aspect_rad <- aspect * (pi / 180)
northness  <- cos(aspect_rad);  names(northness) <- "Northness"
eastness   <- sin(aspect_rad);  names(eastness)  <- "Eastness"

chm       <- rast(raster_path("chm.tif"));        names(chm)       <- "CHM"
dist_cost <- rast(raster_path("dist_coast.tif"));  names(dist_cost) <- "Dist_Coast"
dmcls     <- rast(raster_path("dmcls.tif"));       names(dmcls)     <- "DMCLS"
dnst      <- rast(raster_path("dnst.tif"));        names(dnst)      <- "DNST"

continuous_stack <- c(dem, slope, northness, eastness, chm, dist_cost, dmcls, dnst)
cat(sprintf("  Continuous stack: %d layers\n", nlyr(continuous_stack)))

# TPI reclassification (Weiss, 2001)
tpi_300_cls  <- reclassify_tpi(rast(raster_path("tpi_r300.tif")),  slope, "TPI300_CLS")
tpi_1000_cls <- reclassify_tpi(rast(raster_path("tpi_r1000.tif")), slope, "TPI1000_CLS")
tpi_3000_cls <- reclassify_tpi(rast(raster_path("tpi_r3000.tif")), slope, "TPI3000_CLS")

frtp <- rast(raster_path("FRTP.tif")); names(frtp) <- "FRTP"
categorical_stack <- c(frtp, tpi_300_cls, tpi_1000_cls, tpi_3000_cls)
cat(sprintf("  Categorical stack: %d layers\n", nlyr(categorical_stack)))

# =============================================================================
# 2. Extract per season
# =============================================================================
cat("\n[2/4] Extracting zonal statistics per season...\n")

season_files <- c(
  Spring = "voronoi_aniso_spring.gpkg",
  Summer = "voronoi_aniso_summer.gpkg",
  Autumn = "voronoi_aniso_autumn.gpkg",
  Winter = "voronoi_aniso_winter.gpkg"
)

all_seasons <- list()

for (s in names(season_files)) {
  cat(sprintf("  %s: ", s))
  
  voronoi_sf <- st_read(vector_path(season_files[s]), quiet = TRUE) %>%
    mutate(site = as.character(site))
  
  # Continuous: mean, stdev, max, min
  ext_cont <- exact_extract(continuous_stack, voronoi_sf,
                            c("mean", "stdev", "max", "min"), progress = FALSE)
  
  # Categorical: area fraction
  ext_frac <- exact_extract(categorical_stack, voronoi_sf, "frac", progress = FALSE)
  colnames(ext_frac) <- gsub("frac_", "", colnames(ext_frac))
  colnames(ext_frac) <- paste0(colnames(ext_frac), "_frac")
  ext_frac[is.na(ext_frac)] <- 0
  
  # Combine
  season_df <- voronoi_sf %>%
    st_drop_geometry() %>%
    dplyr::select(site) %>%
    bind_cols(ext_cont) %>%
    bind_cols(ext_frac) %>%
    mutate(season = s)
  
  all_seasons[[s]] <- season_df
  cat(sprintf("%d polygons × %d features\n", nrow(season_df), ncol(season_df) - 2))
}

spatial_aniso <- bind_rows(all_seasons)

# =============================================================================
# 3. Directional feature extraction (wind-direction-aware variables)
# =============================================================================
cat("\n[3/5] Extracting directional features...\n")

prev_wind <- readRDS(inter_path("prevailing_wind.rds")) %>%
  mutate(site = as.character(site))

stations_pts <- st_read(vector_path("station_213.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site))

# Directional extraction function
extract_directional <- function(voronoi_sf, stations_sf, prev_wind_df,
                                dem_rast, chm_rast, season_name) {
  results <- list()
  
  for (i in seq_len(nrow(voronoi_sf))) {
    site_id <- as.character(voronoi_sf$site[i])
    poly <- voronoi_sf[i, ]
    
    stn <- stations_sf %>% dplyr::filter(site == site_id)
    if (nrow(stn) == 0) next
    stn_coords <- st_coordinates(stn)
    stn_x <- stn_coords[1, 1]; stn_y <- stn_coords[1, 2]
    
    wd_row <- prev_wind_df %>% dplyr::filter(site == site_id, season == season_name)
    if (nrow(wd_row) == 0) next
    
    upwind_rad <- wd_row$dom_wd[1] * pi / 180
    uw_dx <- -sin(upwind_rad); uw_dy <- -cos(upwind_rad)
    
    tryCatch({
      dem_crop <- terra::mask(terra::crop(dem_rast, terra::vect(poly)), terra::vect(poly))
      chm_crop <- terra::mask(terra::crop(chm_rast, terra::vect(poly)), terra::vect(poly))
      
      dem_cells <- as.data.frame(dem_crop, xy = TRUE, na.rm = TRUE)
      chm_cells <- as.data.frame(chm_crop, xy = TRUE, na.rm = TRUE)
      if (nrow(dem_cells) < 10) next
      
      names(dem_cells)[3] <- "dem_val"; names(chm_cells)[3] <- "chm_val"
      cells <- merge(dem_cells, chm_cells, by = c("x", "y"), all.x = TRUE)
      cells$chm_val[is.na(cells$chm_val)] <- 0
      
      cells$upwind_dist <- (cells$x - stn_x) * uw_dx + (cells$y - stn_y) * uw_dy
      upwind_cells   <- cells[cells$upwind_dist > 0, ]
      downwind_cells <- cells[cells$upwind_dist <= 0, ]
      
      uw_dem <- if (nrow(upwind_cells) > 0) mean(upwind_cells$dem_val, na.rm = TRUE) else NA
      dw_dem <- if (nrow(downwind_cells) > 0) mean(downwind_cells$dem_val, na.rm = TRUE) else NA
      
      grad <- if (nrow(cells) > 10) coef(lm(dem_val ~ upwind_dist, data = cells))[2] else NA
      fetch <- if (nrow(upwind_cells) > 0) max(upwind_cells$upwind_dist) else 0
      uw_chm <- if (nrow(upwind_cells) > 0) mean(upwind_cells$chm_val, na.rm = TRUE) else NA
      uw_max_dem <- if (nrow(upwind_cells) > 0) max(upwind_cells$dem_val, na.rm = TRUE) else NA
      stn_dem <- dem_cells$dem_val[which.min((dem_cells$x - stn_x)^2 + (dem_cells$y - stn_y)^2)]
      
      results[[length(results) + 1]] <- tibble(
        site = site_id, season = season_name,
        dem_updown_diff     = uw_dem - dw_dem,
        along_wind_gradient = as.numeric(grad),
        upwind_fetch        = fetch,
        upwind_chm_mean     = uw_chm,
        barrier_index       = uw_max_dem - stn_dem
      )
    }, error = function(e) {})
  }
  bind_rows(results)
}

dir_list <- list()
for (s in names(season_files)) {
  cat(sprintf("  %s... ", s))
  vor <- st_read(vector_path(season_files[s]), quiet = TRUE) %>%
    mutate(site = as.character(site))
  dir_s <- extract_directional(vor, stations_pts, prev_wind, dem, chm, s)
  dir_list[[s]] <- dir_s
  cat(sprintf("%d stations\n", nrow(dir_s)))
}

dir_features <- bind_rows(dir_list)

# Merge directional features into main dataframe
spatial_aniso <- spatial_aniso %>%
  dplyr::left_join(dir_features, by = c("site", "season"))

cat(sprintf("  Directional features added: %d rows matched\n", sum(!is.na(spatial_aniso$upwind_fetch))))

# =============================================================================
# 4. Post-processing
# =============================================================================
cat("\n[4/5] Post-processing...\n")

# Fix column separator (exact_extract uses '.' but we standardize)
# Keep as-is — downstream scripts handle both formats

# Derived variable: DEM Range
max_col <- grep("^max[._]DEM$", names(spatial_aniso), value = TRUE)[1]
min_col <- grep("^min[._]DEM$", names(spatial_aniso), value = TRUE)[1]

if (!is.na(max_col) && !is.na(min_col)) {
  spatial_aniso$DEM_Range <- spatial_aniso[[max_col]] - spatial_aniso[[min_col]]
  cat(sprintf("  DEM_Range created (%s − %s)\n", max_col, min_col))
} else {
  warning("Could not find max/min DEM columns for DEM_Range")
}

# Fill NA fractions with 0
spatial_aniso <- spatial_aniso %>%
  mutate(across(ends_with("_frac"), ~replace_na(., 0)))

# =============================================================================
# 5. Filter to common stations and save
# =============================================================================
cat("\n[5/5] Filtering to common stations...\n")

spatial_aniso_common <- spatial_aniso %>%
  filter(site %in% COMMON_SITES)

cat(sprintf("  Full:     %d rows (%d stations × %d seasons)\n",
            nrow(spatial_aniso), length(unique(spatial_aniso$site)),
            length(unique(spatial_aniso$season))))
cat(sprintf("  Filtered: %d rows (%d common stations)\n",
            nrow(spatial_aniso_common), length(unique(spatial_aniso_common$site))))
cat(sprintf("  Features: %d columns\n", ncol(spatial_aniso_common)))

# Save both versions
saveRDS(spatial_aniso, inter_path("spatial_features_aniso_all.rds"))
saveRDS(spatial_aniso_common, inter_path("spatial_features_aniso.rds"))

cat(sprintf("\n[Complete]\n"))
cat(sprintf("  → %s (206 stations)\n", inter_path("spatial_features_aniso.rds")))
cat(sprintf("  → %s (213 stations)\n", inter_path("spatial_features_aniso_all.rds")))
