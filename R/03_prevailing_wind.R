##############################################################################
# 03_prevailing_wind.R
#
# Seasonal prevailing wind direction calculation:
#   - Decompose wind speed/direction into u,v components
#   - Vector-average per station × season (Grange, 2014)
#   - Inverse atan2 to recover meteorological direction (0–360°)
#
# Input:  01_intermediate/wind_daily_clean.rds
#         00_data/vector/station_213.gpkg
# Output: 01_intermediate/prevailing_wind.rds
#         01_intermediate/station_seasonal_wind.gpkg
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  03. Seasonal Prevailing Wind Direction\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("[1/3] Loading data...\n")

wind_df    <- readRDS(inter_path("wind_daily_clean.rds"))
station_sf <- st_read(vector_path("station_213.gpkg"), quiet = TRUE)

cat(sprintf("  Wind records: %s\n", format(nrow(wind_df), big.mark = ",")))
cat(sprintf("  Stations: %d\n", nrow(station_sf)))

# =============================================================================
# 2. Vector decomposition and seasonal averaging
# =============================================================================
cat("[2/3] Computing seasonal prevailing wind direction...\n")

wind_df <- wind_df %>%
  mutate(
    site = as.character(site),
    season = month_to_season(month(Date))
  ) %>%
  filter(ws > 0, !is.na(wd)) %>%
  mutate(
    theta_rad = wd * (pi / 180),
    u_comp = -ws * sin(theta_rad),   # East-West component
    v_comp = -ws * cos(theta_rad)    # North-South component
  )

# Station × season vector average
seasonal_wind <- wind_df %>%
  group_by(site, season) %>%
  summarise(
    mean_u  = mean(u_comp, na.rm = TRUE),
    mean_v  = mean(v_comp, na.rm = TRUE),
    mean_ws = mean(ws, na.rm = TRUE),
    n_obs   = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Inverse atan2 → meteorological direction (0–360°)
    dom_wd = (atan2(mean_u, mean_v) * (180 / pi)) + 180,
    dom_wd = ifelse(dom_wd == 0, 360, dom_wd)
  )

cat(sprintf("  Computed: %d station-season combinations\n", nrow(seasonal_wind)))

# Summary
cat("\n  Prevailing wind direction summary:\n")
seasonal_wind %>%
  group_by(season) %>%
  summarise(
    mean_wd = round(mean(dom_wd), 1),
    sd_wd   = round(sd(dom_wd), 1),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# 3. Merge with spatial data and save
# =============================================================================
cat("\n[3/3] Saving...\n")

# Wide format for spatial join
seasonal_wide <- seasonal_wind %>%
  select(site, season, dom_wd) %>%
  tidyr::pivot_wider(names_from = season, values_from = dom_wd, names_prefix = "wd_")

station_seasonal_sf <- station_sf %>%
  mutate(site = as.character(site)) %>%
  left_join(seasonal_wide, by = "site")

# Save
saveRDS(seasonal_wind, inter_path("prevailing_wind.rds"))
st_write(station_seasonal_sf, inter_path("station_seasonal_wind.gpkg"),
         append = FALSE, quiet = TRUE)

cat(sprintf("  → %s\n", inter_path("prevailing_wind.rds")))
cat(sprintf("  → %s\n", inter_path("station_seasonal_wind.gpkg")))
cat("\n[Complete]\n")
