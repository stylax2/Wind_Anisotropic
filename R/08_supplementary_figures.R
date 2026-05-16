##############################################################################
# 08_supplementary_figures.R
#
# Main text:
#   Fig 3:  Anisotropic vs Isotropic Voronoi comparison map
#   Fig 4:  Seasonal wind speed/direction spatial distribution
#
# Supplementary:
#   Fig S1: Distance-band time series (wind speed + direction arrows)
#   Fig S2: Seasonal wind speed histogram (leeward/windward)
#   Fig S3: Seasonal wind rose
#
# Input:  00_data/vector/voronoi_aniso_*.gpkg, voronoi_iso.gpkg
#         00_data/vector/area_boundary.gpkg, mountain_ridge.gpkg
#         01_intermediate/wind_daily_clean.rds
#         00_data/weather/station_meta.csv
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  08. Supplementary Figures & Voronoi Comparison\n")
cat("================================================================\n\n")

# =============================================================================
# Common data loading
# =============================================================================
cat("[0] Loading common data...\n")

wind_daily <- readRDS(inter_path("wind_daily_clean.rds")) %>%
  mutate(site = as.character(site))

boundary <- st_read(vector_path("area_boundary.gpkg"), quiet = TRUE)

# Mountain ridge (optional)
has_ridge <- file.exists(vector_path("mountain_ridge.gpkg"))
if (has_ridge) {
  ridge <- st_read(vector_path("mountain_ridge.gpkg"), quiet = TRUE)
  cat("  Mountain ridge: loaded\n")
}

stations_sf <- st_read(vector_path("station_213.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site)) %>%
  left_join(STATION_META %>% dplyr::select(site, dist_km, dist_band), by = "site")

cat(sprintf("  Daily wind: %s rows\n", format(nrow(wind_daily), big.mark = ",")))

# =============================================================================
# Fig 3: Anisotropic vs Isotropic Voronoi comparison
# =============================================================================
cat("\n[1/5] Fig 3: Voronoi comparison map...\n")

voronoi_aniso <- st_read(vector_path("voronoi_aniso_spring.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site), type = "Anisotropic (Spring)")
voronoi_iso <- st_read(vector_path("voronoi_iso.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site), type = "Isotropic")

fig3a <- ggplot() +
  geom_sf(data = boundary, fill = "gray95", color = "gray60", linewidth = 0.3) +
  geom_sf(data = voronoi_aniso, fill = NA, color = "steelblue", linewidth = 0.25, alpha = 0.8) +
  { if (has_ridge) geom_sf(data = ridge, color = "darkgreen", linewidth = 0.8, linetype = "solid") } +
  geom_sf(data = stations_sf, shape = 16, size = 1.5, color = "red3", alpha = 0.9) +
  labs(title = "(a) Anisotropic cost-distance Voronoi (Spring)") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.title = element_text(face = "bold", size = 10))

fig3b <- ggplot() +
  geom_sf(data = boundary, fill = "gray95", color = "gray60", linewidth = 0.3) +
  geom_sf(data = voronoi_iso, fill = NA, color = "darkorange", linewidth = 0.25, alpha = 0.8) +
  { if (has_ridge) geom_sf(data = ridge, color = "darkgreen", linewidth = 0.8, linetype = "solid") } +
  geom_sf(data = stations_sf, shape = 16, size = 1.5, color = "red3", alpha = 0.9) +
  labs(title = "(b) Isotropic Voronoi (Thiessen polygon)") +
  theme_void(base_family = FONT_FAMILY) +
  theme(plot.title = element_text(face = "bold", size = 10))

fig3 <- fig3a | fig3b

save_figure(fig3, fig_path("Fig3_voronoi_comparison.tiff"), width = 280, height = 140)

# =============================================================================
# Fig 4: Seasonal wind speed/direction map (4 panels)
# =============================================================================
cat("[2/5] Fig 4: Seasonal wind map...\n")

# Load prevailing wind for arrows
prev_wind <- readRDS(inter_path("prevailing_wind.rds")) %>%
  mutate(site = as.character(site))

season_files <- c(Spring = "voronoi_aniso_spring.gpkg",
                  Summer = "voronoi_aniso_summer.gpkg",
                  Autumn = "voronoi_aniso_autumn.gpkg",
                  Winter = "voronoi_aniso_winter.gpkg")

# Wind p95 for coloring
wind_p95_all <- readRDS(inter_path("wind_seasonal_p95_all.rds"))

fig4_panels <- list()

for (s in SEASONS) {
  vor <- st_read(vector_path(season_files[s]), quiet = TRUE) %>%
    mutate(site = as.character(site)) %>%
    left_join(wind_p95_all %>% filter(season == s) %>% dplyr::select(site, p95_ws),
              by = "site")
  
  # Station centroids with wind arrows
  stn_wind <- stations_sf %>%
    inner_join(prev_wind %>% filter(season == s) %>% dplyr::select(site, dom_wd, mean_ws),
               by = "site") %>%
    mutate(
      lon = st_coordinates(.)[, 1],
      lat = st_coordinates(.)[, 2],
      arrow_len = 2500,
      dx = -sin(dom_wd * pi / 180) * arrow_len,
      dy = -cos(dom_wd * pi / 180) * arrow_len
    )
  
  p <- ggplot() +
    geom_sf(data = vor, aes(fill = p95_ws), color = "gray40", linewidth = 0.1) +
    scale_fill_viridis_c(option = "inferno", name = expression("p95 (m s"^{-1}*")"),
                         limits = c(0, NA), na.value = "gray90") +
    { if (has_ridge) geom_sf(data = ridge, color = "green3", linewidth = 0.6) } +
    geom_segment(data = st_drop_geometry(stn_wind),
                 aes(x = lon, y = lat, xend = lon + dx, yend = lat + dy),
                 arrow = arrow(length = unit(1, "mm"), type = "closed"),
                 color = "dodgerblue", linewidth = 0.3, alpha = 0.7) +
    labs(title = s) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = if (s == "Winter") "right" else "none",
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(0.3, "cm")
    )
  
  fig4_panels[[s]] <- p
}

fig4 <- (fig4_panels$Spring | fig4_panels$Summer) /
  (fig4_panels$Autumn | fig4_panels$Winter)

save_figure(fig4, fig_path("Fig4_seasonal_wind_map.tiff"), width = 260, height = 260)

# =============================================================================
# Fig S1: Distance-band time series
# =============================================================================
cat("[3/5] Fig S1: Distance-band time series...\n")

# Monthly aggregation
wind_monthly <- wind_daily %>%
  mutate(ym = floor_date(Date, "month")) %>%
  group_by(site, ym) %>%
  summarise(
    p95_ws = quantile(ws, 0.95, na.rm = TRUE),
    dom_wd = {
      valid <- !is.na(wd) & ws > 0
      theta <- wd[valid] * pi / 180
      if (length(theta) == 0) NA_real_
      else (atan2(mean(sin(theta)), mean(cos(theta))) * 180 / pi) %% 360
    },
    .groups = "drop"
  ) %>%
  left_join(STATION_META %>% dplyr::select(site, dist_band, dist_km), by = "site") %>%
  filter(!is.na(dist_band))

# Band-level monthly summary
band_monthly <- wind_monthly %>%
  group_by(dist_band, ym) %>%
  summarise(
    mean_p95 = mean(p95_ws, na.rm = TRUE),
    q25_p95  = quantile(p95_ws, 0.25, na.rm = TRUE),
    q75_p95  = quantile(p95_ws, 0.75, na.rm = TRUE),
    mean_wd  = {
      valid <- !is.na(dom_wd)
      theta <- dom_wd[valid] * pi / 180
      if (length(theta) == 0) NA_real_
      else (atan2(mean(sin(theta)), mean(cos(theta))) * 180 / pi) %% 360
    },
    n_stations = n(),
    .groups = "drop"
  )

# Band labels with station counts
band_counts <- STATION_META %>%
  filter(!is.na(dist_band)) %>%
  group_by(dist_band) %>% tally()
band_labels <- setNames(
  paste0(band_counts$dist_band, " (n=", band_counts$n, ")"),
  band_counts$dist_band
)

# Season background shading
make_season_bg <- function(years, season_name, months) {
  do.call(rbind, lapply(years, function(y) {
    data.frame(
      xmin = as.Date(sprintf("%d-%02d-01", y, months[1])),
      xmax = as.Date(sprintf("%d-%02d-28", y, tail(months, 1))),
      season = season_name
    )
  }))
}

season_bg <- rbind(
  make_season_bg(2021:2025, "Winter", c(12, 1, 2)),
  make_season_bg(2021:2025, "Summer", c(6, 7, 8))
)

# (a) Wind speed
pS1a <- ggplot(band_monthly, aes(x = ym, y = mean_p95,
                                 color = dist_band, fill = dist_band)) +
  geom_rect(data = season_bg %>% filter(season == "Winter"),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#DCEEFB", alpha = 0.4) +
  geom_rect(data = season_bg %>% filter(season == "Summer"),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#FFF8DC", alpha = 0.4) +
  geom_ribbon(aes(ymin = q25_p95, ymax = q75_p95), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = BAND_COLORS, labels = band_labels, name = "Distance to coast") +
  scale_fill_manual(values = BAND_COLORS, labels = band_labels, name = "Distance to coast") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.02, 0)) +
  labs(title = "(a) Monthly extreme wind speed by distance-to-coast band",
       subtitle = "Line = band mean; ribbon = IQR (Q25\u2013Q75); blue/yellow shading = winter/summer",
       x = "", y = expression("p95 wind speed (m s"^{-1}*")"))

# (b) Wind direction arrows
band_levels <- c("0\u201320 km", "20\u201340 km", "40\u201360 km", "60\u201380 km", "80+ km")

band_monthly_arrow <- band_monthly %>%
  mutate(
    y_num = as.numeric(factor(dist_band, levels = band_levels)),
    x_num = as.numeric(ym),
    wd_rad = mean_wd * pi / 180,
    arrow_len_x = 10, arrow_len_y = 0.25,
    dx = -sin(wd_rad) * arrow_len_x,
    dy = -cos(wd_rad) * arrow_len_y,
    x_end = x_num + dx,
    y_end = y_num + dy
  )

season_bg_num <- season_bg %>%
  mutate(xmin_n = as.numeric(xmin), xmax_n = as.numeric(xmax))

pS1b <- ggplot(band_monthly_arrow) +
  geom_rect(data = season_bg_num %>% filter(season == "Winter"),
            aes(xmin = xmin_n, xmax = xmax_n, ymin = 0.4, ymax = 5.6),
            fill = "#DCEEFB", alpha = 0.4) +
  geom_rect(data = season_bg_num %>% filter(season == "Summer"),
            aes(xmin = xmin_n, xmax = xmax_n, ymin = 0.4, ymax = 5.6),
            fill = "#FFF8DC", alpha = 0.4) +
  geom_segment(aes(x = x_num, y = y_num, xend = x_end, yend = y_end, color = mean_p95),
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.5, alpha = 0.85) +
  scale_color_viridis_c(option = "inferno", name = expression("p95 (m s"^{-1}*")")) +
  scale_x_continuous(
    breaks = as.numeric(as.Date(paste0(2021:2026, "-01-01"))),
    labels = 2021:2026, expand = c(0.02, 0)) +
  scale_y_continuous(breaks = 1:5, labels = rev(band_levels),
                     limits = c(0.3, 5.7), expand = c(0, 0)) +
  labs(title = "(b) Monthly prevailing wind direction by distance-to-coast band",
       subtitle = "Arrows = direction wind blows toward (vector-averaged); color = band-mean p95",
       x = "Date", y = "Distance to coast")

figS1 <- pS1a / pS1b + plot_layout(heights = c(2, 1))

save_figure(figS1, sup_path("FigS1_distance_band_timeseries.tiff"), width = 280, height = 200)

# =============================================================================
# Fig S2: Seasonal wind speed histogram (leeward/windward)
# =============================================================================
cat("[4/5] Fig S2: Wind speed histogram...\n")

wind_hist <- wind_daily %>%
  left_join(STATION_META %>% dplyr::select(site, dist_km), by = "site") %>%
  mutate(
    season = month_to_season(month(Date)),
    season = factor(season, levels = SEASONS),
    region = ifelse(dist_km <= 45, "Leeward (coastal)", "Windward (inland)")
  ) %>%
  filter(!is.na(region))

# p95 lines
p95_lines <- wind_hist %>%
  group_by(season, region) %>%
  summarise(p95 = quantile(ws, 0.95, na.rm = TRUE), .groups = "drop")

figS2 <- ggplot(wind_hist, aes(x = ws, fill = region)) +
  geom_histogram(bins = 50, alpha = 0.55, position = "identity") +
  geom_vline(data = p95_lines, aes(xintercept = p95, color = region),
             linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~season, scales = "free_y") +
  scale_fill_manual(name = "Region",
                    values = c("Leeward (coastal)" = "#D55E00",
                               "Windward (inland)" = "#0072B2")) +
  scale_color_manual(name = "Region",
                     values = c("Leeward (coastal)" = "#D55E00",
                                "Windward (inland)" = "#0072B2")) +
  labs(x = expression("Daily maximum wind speed (m s"^{-1}*")"),
       y = "Frequency",
       title = "Figure S2. Seasonal wind speed distribution",
       subtitle = "Dashed lines = 95th percentile by region")

save_figure(figS2, sup_path("FigS2_wind_histogram.tiff"), width = 250, height = 180)

# =============================================================================
# Fig S3: Seasonal wind rose
# =============================================================================
cat("[5/5] Fig S3: Wind rose...\n")

wind_rose_data <- wind_daily %>%
  left_join(STATION_META %>% dplyr::select(site, dist_km), by = "site") %>%
  mutate(
    season = month_to_season(month(Date)),
    region = ifelse(dist_km <= 45, "Leeward", "Windward")
  ) %>%
  rename(date = Date) %>%
  filter(!is.na(ws), !is.na(wd), ws > 0)

# Use openair::windRose
png_path <- sup_path("FigS3_wind_rose.png")
tiff_path <- sup_path("FigS3_wind_rose.tiff")

tiff(tiff_path, width = 250, height = 200, units = "mm", res = 600, compression = "lzw",
     family = FONT_FAMILY)
windRose(wind_rose_data, type = "season",
         ws = "ws", wd = "wd",
         cols = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027"),
         key.position = "right",
         paddle = FALSE,
         main = "Figure S3. Seasonal wind rose (all 213 stations)")
dev.off()

# PNG preview
png(png_path, width = 250, height = 200, units = "mm", res = 300, family = FONT_FAMILY)
windRose(wind_rose_data, type = "season",
         ws = "ws", wd = "wd",
         cols = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027"),
         key.position = "right",
         paddle = FALSE,
         main = "Figure S3. Seasonal wind rose (all 213 stations)")
dev.off()

cat(sprintf("  \u2192 %s (TIFF + PNG)\n", basename(tiff_path)))

# =============================================================================
# Summary
# =============================================================================
cat("\n================================================================\n")
cat("  08 Complete\n")
cat("================================================================\n\n")

cat("[Main text]\n")
cat(sprintf("  Fig 3:  %s\n", fig_path("Fig3_voronoi_comparison.tiff")))
cat(sprintf("  Fig 4:  %s\n", fig_path("Fig4_seasonal_wind_map.tiff")))

cat("\n[Supplementary]\n")
cat(sprintf("  Fig S1: %s\n", sup_path("FigS1_distance_band_timeseries.tiff")))
cat(sprintf("  Fig S2: %s\n", sup_path("FigS2_wind_histogram.tiff")))
cat(sprintf("  Fig S3: %s\n", sup_path("FigS3_wind_rose.tiff")))