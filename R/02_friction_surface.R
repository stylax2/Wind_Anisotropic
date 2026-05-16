##############################################################################
# 02_friction_surface.R
#
# Seasonal friction surface generation:
#   - Multi-scale TPI (300m, 1000m, 3000m) → topographic cost
#   - CHM × FRTP seasonal weights → vegetation roughness
#   - Final friction = topographic cost × vegetation roughness
#
# Input:  00_data/raster/ (dem, slope, chm, tpi_r*, FRTP)
# Output: 01_intermediate/friction_spring.tif ... friction_winter.tif
#
# References:
#   Winstral et al. (2002) - multi-scale TPI
#   Belcher et al. (2012)  - slope flow separation
#   Wagenbrenner et al. (2016) - elevation ABL effect
#   Raupach (1994) - z0 ≈ 0.1 × CHM
#   Dolman (1986), Nakai et al. (2008) - seasonal leaf-off roughness
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  02. Seasonal Friction Surface Generation\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load rasters
# =============================================================================
cat("[1/4] Loading rasters...\n")

chm      <- rast(raster_path("chm.tif"))
frtp     <- rast(raster_path("FRTP.tif"))
dem      <- rast(raster_path("dem.tif"))
slope    <- rast(raster_path("slope.tif"))
tpi_300  <- rast(raster_path("tpi_r300.tif"))
tpi_1000 <- rast(raster_path("tpi_r1000.tif"))
tpi_3000 <- rast(raster_path("tpi_r3000.tif"))

# =============================================================================
# 2. Topographic cost (season-invariant)
# =============================================================================
cat("[2/4] Computing topographic cost...\n")

# Z-score normalization
scale_raster <- function(r) {
  r_mean <- global(r, "mean", na.rm = TRUE)[[1]]
  r_sd   <- global(r, "sd",   na.rm = TRUE)[[1]]
  (r - r_mean) / r_sd
}

# Multi-scale TPI combination (Winstral et al., 2002)
tpi_combined_z <- 0.2 * scale_raster(tpi_300) +
                  0.6 * scale_raster(tpi_1000) +
                  0.2 * scale_raster(tpi_3000)

tpi_cost <- exp(-0.2 * tpi_combined_z)

# Slope effect (Belcher et al., 2012)
slope_rad  <- slope * (pi / 180)
slope_cost <- 1 + 2 * tan(slope_rad)

# Elevation effect (Wagenbrenner et al., 2016)
dem_max  <- global(dem, "max", na.rm = TRUE)[[1]]
dem_cost <- 1 - 0.2 * (dem / dem_max)

# Combined topographic cost
topo_cost <- tpi_cost * slope_cost * dem_cost

# =============================================================================
# 3. Seasonal vegetation roughness function
# =============================================================================
cat("[3/4] Building seasonal friction surfaces...\n")

create_friction <- function(chm, frtp_layer, topo_cost, season) {
  base_z0 <- chm * 0.1  # Raupach (1994)

  # FRTP seasonal weights (Dolman 1986; Nakai et al. 2008)
  if (season == "Summer") {
    rcl <- matrix(c(
      1,1.0, 4,1.0, 7,1.0, 10,1.0, 11,1.0, 13,1.0,  # Evergreen conifer
      2,1.0, 8,1.0, 9,1.0,                            # Deciduous broadleaf
      5,1.0, 6,1.0,                                   # Deciduous conifer
      3,1.0,                                          # Mixed
      0,0.1, 12,0.1, 14,0.1                           # Non-forest
    ), ncol = 2, byrow = TRUE)
  } else if (season == "Winter") {
    rcl <- matrix(c(
      1,1.0, 4,1.0, 7,1.0, 10,1.0, 11,1.0, 13,1.0,
      2,0.5, 8,0.5, 9,0.5,                            # Leaf-off: 50% reduction
      5,0.6, 6,0.6,
      3,0.75,
      0,0.1, 12,0.1, 14,0.1
    ), ncol = 2, byrow = TRUE)
  } else {  # Spring & Autumn
    rcl <- matrix(c(
      1,1.0, 4,1.0, 7,1.0, 10,1.0, 11,1.0, 13,1.0,
      2,0.8, 8,0.8, 9,0.8,                            # Transitional
      5,0.8, 6,0.8,
      3,0.9,
      0,0.1, 12,0.1, 14,0.1
    ), ncol = 2, byrow = TRUE)
  }

  weight_layer <- classify(frtp_layer, rcl)
  veg_z0 <- base_z0 * weight_layer
  final_cost <- veg_z0 * topo_cost
  final_cost <- ifel(final_cost <= 0.01, 0.01, final_cost)
  return(final_cost)
}

# =============================================================================
# 4. Generate and save 4 seasonal friction surfaces
# =============================================================================
for (s in SEASONS) {
  cat(sprintf("  %s...", s))
  friction <- create_friction(chm, frtp, topo_cost, s)
  out_file <- inter_path(sprintf("friction_%s.tif", tolower(s)))
  writeRaster(friction, out_file, overwrite = TRUE)
  cat(sprintf(" → %s\n", basename(out_file)))
}

cat("\n[Complete] 4 seasonal friction surfaces saved.\n")
