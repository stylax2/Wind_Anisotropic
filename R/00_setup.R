##############################################################################
# 00_setup.R
#
# Common configuration for the wind-landscape analysis pipeline.
# This script defines:
#   - Required packages (with auto-install)
#   - Directory paths
#   - Common ggplot2 theme and color palettes
#   - 206 common station filtering logic
#   - Utility functions
#
# Usage: source("R/00_setup.R") at the top of every script
##############################################################################

cat("================================================================\n")
cat("  Wind-Landscape Interaction Analysis — Setup\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Package management
# =============================================================================

required_packages <- c(
  # Spatial
  "terra", "sf", "exactextractr", "gdistance",
  # Data manipulation
  "dplyr", "tidyr", "lubridate", "readr",
  # Machine learning
  "xgboost", "ranger", "caret",
  # SHAP
  "SHAPforxgboost",
  # Visualization
  "ggplot2", "patchwork", "ggExtra", "scales",
  # Font
  "extrafont",
  # Wind rose (supplementary)
  "openair"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    cat(sprintf("  Installing: %s\n", paste(missing, collapse = ", ")))
    install.packages(missing, quiet = TRUE)
  }
}

cat("[Packages]\n")
install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))
cat("  All packages loaded.\n")

# =============================================================================
# 2. Directory paths
# =============================================================================

cat("\n[Directories]\n")

# Auto-detect project root (works if sourced from R/ or project root)
if (file.exists("00_data")) {
  PROJECT_ROOT <- getwd()
} else if (file.exists("../00_data")) {
  PROJECT_ROOT <- normalizePath("..")
} else {
  stop("Cannot find 00_data/. Set working directory to project root.")
}

DIR_DATA_RASTER  <- file.path(PROJECT_ROOT, "00_data", "raster")
DIR_DATA_VECTOR  <- file.path(PROJECT_ROOT, "00_data", "vector")
DIR_DATA_WEATHER <- file.path(PROJECT_ROOT, "00_data", "weather")
DIR_INTERMEDIATE <- file.path(PROJECT_ROOT, "01_intermediate")
DIR_FIGURES      <- file.path(PROJECT_ROOT, "02_output", "figures")
DIR_SUPPLEMENT   <- file.path(PROJECT_ROOT, "02_output", "supplementary")

# Create output directories if they don't exist
for (d in c(DIR_INTERMEDIATE, DIR_FIGURES, DIR_SUPPLEMENT)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat(sprintf("  Project root: %s\n", PROJECT_ROOT))

# Verify input data
n_raster <- length(list.files(DIR_DATA_RASTER, pattern = "\\.tif$"))
n_vector <- length(list.files(DIR_DATA_VECTOR, pattern = "\\.gpkg$"))
n_weather <- length(list.files(DIR_DATA_WEATHER))

cat(sprintf("  Rasters: %d files\n", n_raster))
cat(sprintf("  Vectors: %d files\n", n_vector))
cat(sprintf("  Weather: %d files\n", n_weather))

if (n_raster < 11) warning("Expected 11 raster files, found ", n_raster)
if (n_weather < 2) warning("Expected 2 weather files, found ", n_weather)

# =============================================================================
# 3. Path helper functions
# =============================================================================

# Construct full paths from short names
raster_path  <- function(f) file.path(DIR_DATA_RASTER, f)
vector_path  <- function(f) file.path(DIR_DATA_VECTOR, f)
weather_path <- function(f) file.path(DIR_DATA_WEATHER, f)
inter_path   <- function(f) file.path(DIR_INTERMEDIATE, f)
fig_path     <- function(f) file.path(DIR_FIGURES, f)
sup_path     <- function(f) file.path(DIR_SUPPLEMENT, f)

# =============================================================================
# 4. Common constants
# =============================================================================

# Seasons
SEASONS      <- c("Spring", "Summer", "Autumn", "Winter")
SEASON_MONTHS <- list(Spring = 3:5, Summer = 6:8, Autumn = 9:11, Winter = c(12, 1, 2))

# Analysis parameters
TRAIN_RATIO  <- 0.8
RANDOM_SEED  <- 42
N_STATIONS   <- 213
TARGET_VAR   <- "p95_ws"

# =============================================================================
# 5. Font registration and ggplot2 theme
# =============================================================================

cat("\n[Font]\n")

# Ensure Arial is available for PDF/TIFF output
if (requireNamespace("extrafont", quietly = TRUE)) {
  library(extrafont)
  if (!"Arial" %in% fonts()) {
    cat("  Registering system fonts (first run only)...\n")
    font_import(prompt = FALSE)
    loadfonts(device = "win", quiet = TRUE)   # Windows
    loadfonts(device = "pdf", quiet = TRUE)    # PDF
  }
  cat("  Arial loaded via extrafont.\n")
} else if (requireNamespace("showtext", quietly = TRUE)) {
  library(showtext)
  showtext_auto()
  if (!any(grepl("Arial", font_families()))) {
    font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf",
             italic = "ariali.ttf")
  }
  cat("  Arial loaded via showtext.\n")
} else {
  cat("  [Note] Install 'extrafont' or 'showtext' for reliable Arial rendering.\n")
  cat("         install.packages('extrafont')  # recommended\n")
}

FONT_FAMILY <- "Arial"

theme_paper <- theme_bw(base_size = 11, base_family = FONT_FAMILY) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    legend.title  = element_text(face = "bold", size = 9),
    legend.text   = element_text(size = 8),
    strip.background = element_rect(fill = "gray95"),
    strip.text    = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

# Set default theme for all plots
theme_set(theme_paper)

# Set default geom text/label font
update_geom_defaults("text",  list(family = FONT_FAMILY))
update_geom_defaults("label", list(family = FONT_FAMILY))

cat(sprintf("  Font family: %s\n", FONT_FAMILY))
cat("  theme_paper set as default for all ggplot2 output.\n")

# Season colors
SEASON_COLORS <- c(
  Spring = "#2CA02C",
  Summer = "#D62728",
  Autumn = "#FF7F0E",
  Winter = "#1F77B4"
)

# Station type colors
STATION_COLORS <- c(
  ASOS = "#E41A1C",
  AWS  = "#377EB8",
  FMS  = "#4DAF4A",
  NPMA = "#984EA3"
)

# Distance-to-coast band colors (coastal → inland)
BAND_COLORS <- c(
  "0–20 km"  = "#D62728",
  "20–40 km" = "#FF7F0E",
  "40–60 km" = "#9467BD",
  "60–80 km" = "#2CA02C",
  "80+ km"   = "#1F77B4"
)

# SHAP group colors
SHAP_GROUP_COLORS <- list(
  elevation = c("< 600 m a.s.l." = "#4393C3", "≥ 600 m a.s.l." = "#D6604D"),
  phenology = c("Leaf-off (Winter/Spring)" = "#0072B2", "Leaf-on (Summer)" = "#D55E00"),
  aspect    = c("Leeward (East-facing)" = "#D62728", "Other aspects" = "#1F77B4")
)

# =============================================================================
# 6. Common station ID identification (206 common sites)
# =============================================================================

cat("\n[Station filtering]\n")

# Load station metadata
STATION_META <- read.csv(weather_path("station_meta.csv"),
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8") %>%
  mutate(
    site = gsub('"', '', as.character(site)),
    dist_km = dist_coast / 1000,
    dist_band = cut(dist_km,
                    breaks = c(0, 20, 40, 60, 80, Inf),
                    labels = c("0–20 km", "20–40 km", "40–60 km",
                               "60–80 km", "80+ km"),
                    include.lowest = TRUE)
  )

cat(sprintf("  Total stations in metadata: %d\n", nrow(STATION_META)))

# Identify common sites between anisotropic and isotropic Voronoi
# Anisotropic: use spring as reference
if (file.exists(vector_path("voronoi_aniso_spring.gpkg"))) {
  aniso_sites <- st_read(vector_path("voronoi_aniso_spring.gpkg"), quiet = TRUE) %>%
    st_drop_geometry() %>%
    mutate(site = as.character(site)) %>%
    pull(site) %>% unique()
} else {
  aniso_sites <- STATION_META$site
  cat("  [Warning] voronoi_aniso_spring.gpkg not found; using all stations.\n")
}

# Isotropic
if (file.exists(vector_path("voronoi_iso.gpkg"))) {
  iso_sites <- st_read(vector_path("voronoi_iso.gpkg"), quiet = TRUE) %>%
    st_drop_geometry() %>%
    mutate(site = as.character(site)) %>%
    pull(site) %>% unique()
} else {
  iso_sites <- STATION_META$site
  cat("  [Warning] voronoi_iso.gpkg not found; will be created in 05b.\n")
}

COMMON_SITES <- intersect(aniso_sites, iso_sites)
EXCLUDED_SITES <- setdiff(union(aniso_sites, iso_sites), COMMON_SITES)

cat(sprintf("  Anisotropic stations: %d\n", length(aniso_sites)))
cat(sprintf("  Isotropic stations:   %d\n", length(iso_sites)))
cat(sprintf("  Common stations:      %d\n", length(COMMON_SITES)))

if (length(EXCLUDED_SITES) > 0) {
  cat(sprintf("  Excluded stations:    %d (%s)\n",
              length(EXCLUDED_SITES), paste(EXCLUDED_SITES, collapse = ", ")))
}

# Save for downstream scripts
saveRDS(COMMON_SITES, inter_path("common_sites.rds"))

# =============================================================================
# 7. Utility functions
# =============================================================================

# Performance metrics
calc_metrics <- function(obs, pred, label = "Overall") {
  tibble(
    Model = label,
    N     = length(obs),
    RMSE  = round(sqrt(mean((obs - pred)^2)), 3),
    MAE   = round(mean(abs(obs - pred)), 3),
    R2    = round(1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2), 3),
    MAPE  = round(mean(abs((obs - pred) / obs)) * 100, 1)
  )
}

# Seasonal + overall metrics
calc_all_metrics <- function(test_data, pred, label) {
  df <- test_data %>%
    dplyr::select(site, season, p95_ws) %>%
    mutate(pred = pred)
  
  overall <- calc_metrics(df$p95_ws, df$pred, paste0(label, " | All"))
  
  by_season <- df %>%
    group_by(season) %>%
    group_modify(~calc_metrics(.x$p95_ws, .x$pred,
                               paste0(label, " | ", .y$season))) %>%
    ungroup() %>%
    dplyr::select(-season)
  
  bind_rows(overall, by_season)
}

# Weiss (2001) TPI reclassification
reclassify_tpi <- function(tpi_rast, slope_rast, rast_name) {
  t_mean <- global(tpi_rast, "mean", na.rm = TRUE)[1, 1]
  t_sd   <- global(tpi_rast, "sd",   na.rm = TRUE)[1, 1]
  z <- (tpi_rast - t_mean) / t_sd
  
  tpi_cls <- ifel(z <= -1.0, 1,
                  ifel(z > -1.0 & z <= -0.5, 2,
                       ifel(z > -0.5 & z < 0.5 & slope_rast <= 5, 3,
                            ifel(z > -0.5 & z < 0.5 & slope_rast > 5, 4,
                                 ifel(z >= 0.5 & z < 1.0, 5, 6)))))
  names(tpi_cls) <- rast_name
  return(tpi_cls)
}

# Assign month to season
month_to_season <- function(m) {
  case_when(
    m %in% 3:5  ~ "Spring",
    m %in% 6:8  ~ "Summer",
    m %in% 9:11 ~ "Autumn",
    TRUE         ~ "Winter"
  )
}

# Save figure as TIFF (600 dpi) + PNG (300 dpi, preview)
save_figure <- function(plot, filename, width = 180, height = 120, dpi_tiff = 600) {
  # TIFF for journal submission
  ggsave(filename,
         plot = plot,
         width = width, height = height, units = "mm",
         dpi = dpi_tiff, device = "tiff", compression = "lzw")
  
  # PNG for quick preview
  png_file <- sub("\\.tiff?$", ".png", filename)
  ggsave(png_file,
         plot = plot,
         width = width, height = height, units = "mm",
         dpi = 300)
  
  cat(sprintf("  → %s (TIFF %ddpi + PNG)\n", basename(filename), dpi_tiff))
}

# =============================================================================
# 8. Summary
# =============================================================================

cat("\n================================================================\n")
cat("  Setup complete.\n")
cat(sprintf("  Common stations: %d\n", length(COMMON_SITES)))
cat(sprintf("  Target variable: %s\n", TARGET_VAR))
cat(sprintf("  Seasons: %s\n", paste(SEASONS, collapse = ", ")))
cat(sprintf("  Random seed: %d\n", RANDOM_SEED))
cat("================================================================\n\n")
