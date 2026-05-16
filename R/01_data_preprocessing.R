##############################################################################
# 01_data_preprocessing.R
#
# Wind data preprocessing:
#   1. Load daily wind speed/direction data (213 stations × 5 years)
#   2. Quality control: remove missing values, flag anomalies
#   3. Aggregate to station × season level: p95, mean, max wind speed
#   4. Filter to 206 common stations
#
# Input:  00_data/weather/wind_daily.rds
# Output: 01_intermediate/wind_seasonal_p95.rds
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  01. Data Preprocessing\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load raw data
# =============================================================================
cat("[1/4] Loading daily wind data...\n")

wind_raw <- readRDS(weather_path("wind_daily.rds"))

cat(sprintf("  Rows: %s | Columns: %s\n", 
            format(nrow(wind_raw), big.mark = ","), ncol(wind_raw)))
cat(sprintf("  Date range: %s to %s\n", min(wind_raw$Date), max(wind_raw$Date)))
cat(sprintf("  Stations: %d\n", length(unique(wind_raw$site))))

# =============================================================================
# 2. Quality control
# =============================================================================
cat("\n[2/4] Quality control...\n")

wind_clean <- wind_raw %>%
  mutate(site = as.character(site)) %>%
  filter(
    !is.na(ws),           # Remove missing wind speed
    !is.na(wd),           # Remove missing wind direction
    ws >= 0,              # Remove negative values
    wd >= 0 & wd <= 360   # Valid direction range
  )

n_removed <- nrow(wind_raw) - nrow(wind_clean)
cat(sprintf("  Removed: %s rows (%.1f%%)\n",
            format(n_removed, big.mark = ","),
            n_removed / nrow(wind_raw) * 100))
cat(sprintf("  Retained: %s rows\n", format(nrow(wind_clean), big.mark = ",")))

# Per-station observation count
station_obs <- wind_clean %>%
  group_by(site) %>%
  summarise(
    n_days = n(),
    date_min = min(Date),
    date_max = max(Date),
    .groups = "drop"
  )

cat(sprintf("  Observation days per station: %d–%d (median: %d)\n",
            min(station_obs$n_days), max(station_obs$n_days),
            median(station_obs$n_days)))

# =============================================================================
# 3. Seasonal aggregation
# =============================================================================
cat("\n[3/4] Seasonal aggregation (p95 wind speed)...\n")

wind_seasonal <- wind_clean %>%
  mutate(season = month_to_season(month(Date))) %>%
  group_by(site, season) %>%
  summarise(
    p95_ws  = quantile(ws, 0.95, na.rm = TRUE),
    mean_ws = mean(ws, na.rm = TRUE),
    max_ws  = max(ws, na.rm = TRUE),
    sd_ws   = sd(ws, na.rm = TRUE),
    n_obs   = n(),
    .groups = "drop"
  )

cat(sprintf("  Aggregated: %d rows (stations × seasons)\n", nrow(wind_seasonal)))
cat(sprintf("  Stations with all 4 seasons: %d\n",
            wind_seasonal %>% group_by(site) %>% tally() %>% filter(n == 4) %>% nrow()))

# Summary statistics
cat("\n  p95 wind speed summary by season:\n")
p95_summary <- wind_seasonal %>%
  group_by(season) %>%
  summarise(
    mean  = round(mean(p95_ws), 2),
    sd    = round(sd(p95_ws), 2),
    min   = round(min(p95_ws), 2),
    max   = round(max(p95_ws), 2),
    .groups = "drop"
  )
print(p95_summary)

# =============================================================================
# 4. Filter to common stations (206)
# =============================================================================
cat("\n[4/4] Filtering to common stations...\n")

wind_seasonal_common <- wind_seasonal %>%
  filter(site %in% COMMON_SITES)

cat(sprintf("  Before: %d rows (%d stations)\n",
            nrow(wind_seasonal), length(unique(wind_seasonal$site))))
cat(sprintf("  After:  %d rows (%d stations)\n",
            nrow(wind_seasonal_common), length(unique(wind_seasonal_common$site))))

# =============================================================================
# 5. Save
# =============================================================================

# Full 213-station version (for Fig S1, descriptive stats)
saveRDS(wind_seasonal, inter_path("wind_seasonal_p95_all.rds"))

# 206-station version (for modeling)
saveRDS(wind_seasonal_common, inter_path("wind_seasonal_p95.rds"))

# Clean daily data (for supplementary figures)
saveRDS(wind_clean, inter_path("wind_daily_clean.rds"))

cat("\n[Complete]\n")
cat(sprintf("  → %s (206 stations, modeling)\n", inter_path("wind_seasonal_p95.rds")))
cat(sprintf("  → %s (213 stations, descriptive)\n", inter_path("wind_seasonal_p95_all.rds")))
cat(sprintf("  → %s (daily, supplementary figures)\n", inter_path("wind_daily_clean.rds")))
