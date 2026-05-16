##############################################################################
# 05b_zonal_stats_isotropic.R
#
# Isotropic (Euclidean distance) Voronoi comparison:
#   1. Generate isotropic Voronoi from R (or load from QGIS)
#   2. Extract identical zonal statistics as 05a
#   3. Identify common stations (aniso ∩ iso)
#   4. Ensure column-level alignment with anisotropic data
#
# Input:  00_data/vector/station_213.gpkg
#         00_data/vector/area_boundary.gpkg
#         00_data/vector/voronoi_iso.gpkg (optional, from QGIS)
#         00_data/raster/*.tif
# Output: 01_intermediate/spatial_features_iso.rds
#         01_intermediate/common_sites.rds (updated)
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  05b. Zonal Statistics — Isotropic Voronoi\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load or generate isotropic Voronoi
# =============================================================================
cat("[1/5] Isotropic Voronoi polygons...\n")

iso_path <- vector_path("voronoi_iso.gpkg")

if (file.exists(iso_path)) {
  cat("  Loading existing file from QGIS...\n")
  voronoi_iso <- st_read(iso_path, quiet = TRUE) %>%
    mutate(site = as.character(site))
  cat(sprintf("  Loaded: %d polygons\n", nrow(voronoi_iso)))
  
} else {
  cat("  Generating from R (st_voronoi)...\n")
  
  stations <- st_read(vector_path("station_213.gpkg"), quiet = TRUE)
  boundary <- st_read(vector_path("area_boundary.gpkg"), quiet = TRUE)
  
  # Ensure same CRS
  boundary <- st_transform(boundary, st_crs(stations))
  
  # Generate Voronoi
  station_union <- st_union(stations)
  voronoi_raw <- st_voronoi(station_union, envelope = st_as_sfc(st_bbox(boundary)))
  voronoi_sf  <- st_collection_extract(voronoi_raw, "POLYGON") %>%
    st_as_sf() %>%
    st_intersection(st_union(boundary))
  
  # Assign station IDs via spatial join
  voronoi_iso <- st_join(voronoi_sf, stations) %>%
    filter(!is.na(site)) %>%
    mutate(site = as.character(site))
  
  # Save
  st_write(voronoi_iso, iso_path, delete_layer = TRUE, quiet = TRUE)
  cat(sprintf("  Generated and saved: %d polygons\n", nrow(voronoi_iso)))
}

# =============================================================================
# 2. Update common sites
# =============================================================================
cat("\n[2/5] Updating common stations...\n")

iso_sites   <- unique(voronoi_iso$site)
aniso_sites <- unique(
  st_read(vector_path("voronoi_aniso_spring.gpkg"), quiet = TRUE) %>%
    st_drop_geometry() %>%
    mutate(site = as.character(site)) %>%
    pull(site)
)

COMMON_SITES   <- intersect(aniso_sites, iso_sites)
EXCLUDED_SITES <- setdiff(union(aniso_sites, iso_sites), COMMON_SITES)

cat(sprintf("  Anisotropic: %d | Isotropic: %d | Common: %d | Excluded: %d\n",
            length(aniso_sites), length(iso_sites),
            length(COMMON_SITES), length(EXCLUDED_SITES)))

saveRDS(COMMON_SITES, inter_path("common_sites.rds"))

# Filter isotropic polygons to common sites
voronoi_iso <- voronoi_iso %>% filter(site %in% COMMON_SITES)

# =============================================================================
# 3. Load rasters (same as 05a)
# =============================================================================
cat("\n[3/5] Loading rasters...\n")

dem   <- rast(raster_path("dem.tif"));        names(dem)   <- "DEM"
slope <- rast(raster_path("slope.tif"));      names(slope) <- "Slope"

aspect_rad <- rast(raster_path("aspect.tif")) * (pi / 180)
northness  <- cos(aspect_rad);  names(northness) <- "Northness"
eastness   <- sin(aspect_rad);  names(eastness)  <- "Eastness"

chm       <- rast(raster_path("chm.tif"));        names(chm)       <- "CHM"
dist_cost <- rast(raster_path("dist_coast.tif"));  names(dist_cost) <- "Dist_Coast"
dmcls     <- rast(raster_path("dmcls.tif"));       names(dmcls)     <- "DMCLS"
dnst      <- rast(raster_path("dnst.tif"));        names(dnst)      <- "DNST"

continuous_stack <- c(dem, slope, northness, eastness, chm, dist_cost, dmcls, dnst)

tpi_300_cls  <- reclassify_tpi(rast(raster_path("tpi_r300.tif")),  slope, "TPI300_CLS")
tpi_1000_cls <- reclassify_tpi(rast(raster_path("tpi_r1000.tif")), slope, "TPI1000_CLS")
tpi_3000_cls <- reclassify_tpi(rast(raster_path("tpi_r3000.tif")), slope, "TPI3000_CLS")

frtp <- rast(raster_path("FRTP.tif")); names(frtp) <- "FRTP"
categorical_stack <- c(frtp, tpi_300_cls, tpi_1000_cls, tpi_3000_cls)

# =============================================================================
# 4. Extract zonal statistics (1 set, replicated to 4 seasons)
# =============================================================================
cat("\n[4/5] Extracting zonal statistics...\n")

# Isotropic polygons are season-invariant → extract once, replicate 4×
cat("  Continuous variables...\n")
ext_cont <- exact_extract(continuous_stack, voronoi_iso,
                          c("mean", "stdev", "max", "min"), progress = FALSE)

cat("  Categorical variables...\n")
ext_frac <- exact_extract(categorical_stack, voronoi_iso, "frac", progress = FALSE)
colnames(ext_frac) <- gsub("frac_", "", colnames(ext_frac))
colnames(ext_frac) <- paste0(colnames(ext_frac), "_frac")
ext_frac[is.na(ext_frac)] <- 0

# Base dataframe (single extraction)
base_df <- voronoi_iso %>%
  st_drop_geometry() %>%
  dplyr::select(site) %>%
  bind_cols(ext_cont) %>%
  bind_cols(ext_frac)

cat(sprintf("  Extracted: %d stations × %d features\n", nrow(base_df), ncol(base_df) - 1))

# Replicate to 4 seasons
spatial_iso <- bind_rows(
  lapply(SEASONS, function(s) base_df %>% mutate(season = s))
)

# Fill NA fractions
spatial_iso <- spatial_iso %>%
  mutate(across(ends_with("_frac"), ~replace_na(., 0)))

# DEM Range
max_col <- grep("^max[._]DEM$", names(spatial_iso), value = TRUE)[1]
min_col <- grep("^min[._]DEM$", names(spatial_iso), value = TRUE)[1]

if (!is.na(max_col) && !is.na(min_col)) {
  spatial_iso$DEM_Range <- spatial_iso[[max_col]] - spatial_iso[[min_col]]
  cat(sprintf("  DEM_Range created (%s − %s)\n", max_col, min_col))
}

# =============================================================================
# 5. Directional feature extraction (same logic, isotropic polygons)
# =============================================================================
cat("\n[5/7] Extracting directional features from isotropic polygons...\n")

prev_wind <- readRDS(inter_path("prevailing_wind.rds")) %>%
  mutate(site = as.character(site))

stations_pts <- st_read(vector_path("station_213.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site))

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
        dem_updown_diff = uw_dem - dw_dem,
        along_wind_gradient = as.numeric(grad),
        upwind_fetch = fetch,
        upwind_chm_mean = uw_chm,
        barrier_index = uw_max_dem - stn_dem
      )
    }, error = function(e) {})
  }
  bind_rows(results)
}

# Isotropic: same polygon for all seasons, but wind direction changes
dir_iso_list <- list()
for (s in SEASONS) {
  cat(sprintf("  %s... ", s))
  dir_s <- extract_directional(voronoi_iso, stations_pts, prev_wind, dem, chm, s)
  dir_iso_list[[s]] <- dir_s
  cat(sprintf("%d stations\n", nrow(dir_s)))
}

dir_iso <- bind_rows(dir_iso_list)

spatial_iso <- spatial_iso %>%
  dplyr::left_join(dir_iso, by = c("site", "season"))

cat(sprintf("  Directional features added: %d rows matched\n", sum(!is.na(spatial_iso$upwind_fetch))))

# =============================================================================
# 6. Column alignment with anisotropic data
# =============================================================================
cat("\n[6/7] Aligning columns with anisotropic data...\n")

spatial_aniso <- readRDS(inter_path("spatial_features_aniso.rds"))

cols_aniso <- sort(names(spatial_aniso))
cols_iso   <- sort(names(spatial_iso))

missing_in_iso <- setdiff(cols_aniso, cols_iso)
extra_in_iso   <- setdiff(cols_iso, cols_aniso)

if (length(missing_in_iso) > 0) {
  cat(sprintf("  Adding %d missing columns (filled with 0): %s\n",
              length(missing_in_iso),
              paste(head(missing_in_iso, 5), collapse = ", ")))
  for (col in missing_in_iso) spatial_iso[[col]] <- 0
}

if (length(extra_in_iso) > 0) {
  cat(sprintf("  Removing %d extra columns: %s\n",
              length(extra_in_iso),
              paste(head(extra_in_iso, 5), collapse = ", ")))
}

# Match column order
spatial_iso <- spatial_iso %>% dplyr::select(all_of(names(spatial_aniso)))

cat(sprintf("  Anisotropic: %d rows × %d cols\n", nrow(spatial_aniso), ncol(spatial_aniso)))
cat(sprintf("  Isotropic:   %d rows × %d cols\n", nrow(spatial_iso), ncol(spatial_iso)))

# =============================================================================
# Save
# =============================================================================
saveRDS(spatial_iso, inter_path("spatial_features_iso.rds"))

cat(sprintf("\n[Complete]\n"))
cat(sprintf("  → %s\n", inter_path("spatial_features_iso.rds")))
cat(sprintf("  → %s (updated)\n", inter_path("common_sites.rds")))
