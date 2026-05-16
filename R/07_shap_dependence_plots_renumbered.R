##############################################################################
# 07_shap_dependence_plots.R
#
# SHAP-based variable interpretation (★ FIGURE NUMBERS ALIGNED WITH MANUSCRIPT ★)
#
# Manuscript figure numbering (Landscape Ecology submission):
#   Figure 6:  SHAP summary plot (top-15 variable importance)
#   Figure 7:  Along-wind gradient × coast proximity   (SHAP rank 1)
#   Figure 8:  Downwind-upwind DEM diff × eastness     (SHAP rank 3, directional)
#   Figure 9:  Sub-grid topographic complexity × elevation  (turbulence threshold)
#   Figure 10: Broadleaf fraction × phenology          (leaf-on/leaf-off)
#   Figure 11: Distance to coast × eastness            (foehn amplification zone)
#
# Input:  01_intermediate/shap_long.rds
#         01_intermediate/X_train.rds
# Output: 02_output/figures/Fig06_summary.tiff
#         02_output/figures/Fig07_gradient.tiff
#         02_output/figures/Fig08_updown.tiff
#         02_output/figures/Fig09_terrain.tiff
#         02_output/figures/Fig10_phenology.tiff
#         02_output/figures/Fig11_foehn.tiff
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  07. SHAP Dependence & Summary Plots (manuscript numbering)\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load SHAP data
# =============================================================================
cat("[1/8] Loading SHAP data...\n")

shap_long <- readRDS(inter_path("shap_long.rds"))
X_train   <- readRDS(inter_path("X_train.rds"))

all_vars <- unique(shap_long$variable)
cat(sprintf("  Variables: %d\n", length(all_vars)))

# Identify key column names (exact_extract naming)
find_col <- function(pattern) {
  matched <- grep(pattern, all_vars, value = TRUE)[1]
  if (is.na(matched)) warning(sprintf("Column not found: %s", pattern))
  matched
}

dist_col      <- find_col("mean.*Dist.Coast|mean.*Dist_Coast")
max_dist_col  <- find_col("max.*Dist.Coast|max.*Dist_Coast")
east_col      <- find_col("mean.*Eastness|mean.*eastness")
stdev_dem_col <- find_col("stdev.*DEM|stdev\\.DEM")
mean_dem_col  <- find_col("mean.*DEM|mean\\.DEM")
frtp2_col     <- find_col("FRTP.*2.*frac|2.*FRTP.*frac")
season_sum    <- find_col("season.*Summer|Summer")

# Directional variables
gradient_col  <- find_col("along_wind_gradient")
updown_col    <- find_col("dem_updown_diff")
fetch_col     <- find_col("upwind_fetch")
barrier_col   <- find_col("barrier_index")

cat(sprintf("  Coast dist:  %s\n", dist_col))
cat(sprintf("  Eastness:    %s\n", east_col))
cat(sprintf("  DEM stdev:   %s\n", stdev_dem_col))
cat(sprintf("  DEM mean:    %s\n", mean_dem_col))
cat(sprintf("  Broadleaf:   %s\n", frtp2_col))
cat(sprintf("  Summer:      %s\n", season_sum))
cat(sprintf("  Gradient:    %s\n", gradient_col))
cat(sprintf("  Updown diff: %s\n", updown_col))

# =============================================================================
# 2. Helper: extract SHAP pair data
# =============================================================================
extract_shap_pair <- function(shap_long, x_var, color_var) {
  x_data <- shap_long %>%
    filter(variable == x_var) %>%
    dplyr::select(ID, shap_x = value, raw_x = rfvalue)

  c_data <- shap_long %>%
    filter(variable == color_var) %>%
    dplyr::select(ID, raw_color = rfvalue)

  inner_join(x_data, c_data, by = "ID")
}

# =============================================================================
# Figure 6: SHAP Summary Plot (variable importance, top 15)
# =============================================================================
cat("\n[2/8] Figure 6: SHAP summary plot...\n")

imp_top <- shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 15)

fig06_data <- shap_long %>%
  filter(variable %in% imp_top$variable) %>%
  mutate(variable = factor(variable, levels = rev(imp_top$variable)))

fig06 <- ggplot(fig06_data, aes(x = value, y = variable, color = rfvalue)) +
  geom_jitter(alpha = 0.4, size = 0.8, height = 0.2) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Feature\nvalue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "SHAP value (impact on prediction)", y = "") +
  theme(axis.text.y = element_text(size = 10))

save_figure(fig06, fig_path("Fig06_SHAP_summary.tiff"),
            width = 200, height = 180)

# =============================================================================
# Figure 7: Along-wind gradient × coast proximity (SHAP rank 1)
# =============================================================================
cat("[3/8] Figure 7: Along-wind gradient x coast proximity...\n")

if (!is.na(gradient_col) && !is.na(dist_col)) {
  df07 <- extract_shap_pair(shap_long, gradient_col, dist_col) %>%
    mutate(
      dist_km = raw_color / 1000,
      coast_group = factor(
        ifelse(dist_km <= 40,
               "\u2264 40 km (coastal)",
               "> 40 km (inland)"),
        levels = c("\u2264 40 km (coastal)", "> 40 km (inland)")
      )
    )
  
  fig07 <- ggplot(df07, aes(x = raw_x, y = shap_x)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    # Left shading: negative gradient = upwind higher = foehn descent
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "#FF6B6B", alpha = 0.05) +
    # Right shading: positive gradient = upwind lower = ascending flow
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#6B9FFF", alpha = 0.05) +
    geom_point(aes(color = coast_group, shape = coast_group),
               alpha = 0.5, size = 1.8) +
    geom_smooth(aes(color = coast_group, fill = coast_group),
                method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    scale_color_manual(name = "Coast\nproximity",
                       values = c("\u2264 40 km (coastal)" = "#D62728",
                                  "> 40 km (inland)" = "#1F77B4")) +
    scale_fill_manual(name = "Coast\nproximity",
                      values = c("\u2264 40 km (coastal)" = "#D62728",
                                 "> 40 km (inland)" = "#1F77B4")) +
    scale_shape_manual(name = "Coast\nproximity",
                       values = c("\u2264 40 km (coastal)" = 17,
                                  "> 40 km (inland)" = 16)) +
    annotate("text",
             x = min(df07$raw_x, na.rm = TRUE) * 0.7,
             y = max(df07$shap_x, na.rm = TRUE) * 0.85,
             label = "Negative gradient\n(upwind higher\n= descending flow)",
             size = 2.8, color = "#D62728", fontface = "italic") +
    annotate("text",
             x = max(df07$raw_x, na.rm = TRUE) * 0.7,
             y = min(df07$shap_x, na.rm = TRUE) * 0.85,
             label = "Positive gradient\n(upwind lower\n= ascending flow)",
             size = 2.8, color = "#1F77B4", fontface = "italic") +
    labs(
      x = "Along-wind elevation gradient (m/m)",
      y = "SHAP value (contribution to extreme wind speed)"
    )
  
  fig07_margin <- ggExtra::ggMarginal(fig07, type = "histogram", margins = "x",
                                       size = 8, fill = "gray80", color = "gray50")
  
  save_figure(fig07_margin, fig_path("Fig07_SHAP_gradient.tiff"),
              width = 200, height = 150)
} else {
  cat("  [SKIP] along_wind_gradient not found in SHAP data\n")
}

# =============================================================================
# Figure 8: Downwind-upwind DEM diff × eastness (SHAP rank 3, directional)
# ★ LABELS FIXED to match manuscript §2.6, §3.5, §4.2 ★
#
# Variable definition (manuscript §2.6, corrected):
#   dem_updown_diff = downwind_mean - upwind_mean
#   Negative = upwind higher (foehn-favorable descending flow)
#   Positive = downwind higher (upwind-barrier sheltering / upslope)
# =============================================================================
cat("[4/8] Figure 8: Downwind-upwind DEM diff x eastness (labels fixed)...\n")

if (!is.na(updown_col) && !is.na(east_col)) {
  df08 <- extract_shap_pair(shap_long, updown_col, east_col) %>%
    mutate(
      # =========================================================================
      # [Sign check] If stored rfvalue is actually "upwind - downwind"
      # (opposite of §2.6 definition), uncomment the next line to flip sign.
      # Diagnostic code at the bottom of this script tells you which case applies.
      # -------------------------------------------------------------------------
      # raw_x = -raw_x,
      # =========================================================================
      
      aspect_group = factor(
        ifelse(raw_color > 0.3,
               "Leeward (East-facing)",
               "Other aspects"),
        levels = c("Other aspects", "Leeward (East-facing)")
      )
    )
  
  fig08 <- ggplot(df08, aes(x = raw_x, y = shap_x)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    
    # ---- Background shading (★ swapped ★) ----
    # Left (negative) = upwind higher = foehn descent → red tint
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "#FF6B6B", alpha = 0.05) +
    # Right (positive) = downwind higher = sheltered → blue tint
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#6B9FFF", alpha = 0.05) +
    
    geom_point(aes(color = aspect_group, shape = aspect_group),
               alpha = 0.5, size = 1.8) +
    geom_smooth(aes(color = aspect_group, fill = aspect_group),
                method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    
    scale_color_manual(name = "Aspect group",
                       values = SHAP_GROUP_COLORS$aspect) +
    scale_fill_manual(name = "Aspect group",
                      values = SHAP_GROUP_COLORS$aspect) +
    scale_shape_manual(name = "Aspect group",
                       values = c("Other aspects" = 16,
                                  "Leeward (East-facing)" = 17)) +
    
    # ---- Left label (★ swapped ★) ----
    # Negative = upwind higher = foehn descent (red)
    annotate("text",
             x = min(df08$raw_x, na.rm = TRUE) * 0.6,
             y = max(df08$shap_x, na.rm = TRUE) * 0.85,
             label = "Upwind higher\n(mountain barrier\n= foehn descent)",
             size = 2.8, color = "#D62728", fontface = "italic") +
    
    # ---- Right label (★ swapped ★) ----
    # Positive = downwind higher = sheltered valley (blue)
    annotate("text",
             x = max(df08$raw_x, na.rm = TRUE) * 0.6,
             y = max(df08$shap_x, na.rm = TRUE) * 0.85,
             label = "Downwind higher\n(sheltered valley)",
             size = 2.8, color = "#1F77B4", fontface = "italic") +
    
    # ---- Axis label (★ updated ★) ----
    labs(
      x = "Downwind \u2212 upwind mean elevation (m)",
      y = "SHAP value (contribution to extreme wind speed)"
    )
  
  fig08_margin <- ggExtra::ggMarginal(fig08, type = "histogram", margins = "x",
                                       size = 8, fill = "gray80", color = "gray50")
  
  save_figure(fig08_margin, fig_path("Fig08_SHAP_updown.tiff"),
              width = 200, height = 150)
} else {
  cat("  [SKIP] dem_updown_diff not found in SHAP data\n")
}

# =============================================================================
# Figure 9: Sub-grid topographic complexity × mean elevation
#           (Turbulence threshold response)
# =============================================================================
cat("[5/8] Figure 9: Terrain complexity x elevation...\n")

df09 <- extract_shap_pair(shap_long, stdev_dem_col, mean_dem_col) %>%
  mutate(
    elev_group = factor(
      ifelse(raw_color >= 600,
             "\u2265 600 m a.s.l.",
             "< 600 m a.s.l."),
      levels = c("< 600 m a.s.l.", "\u2265 600 m a.s.l.")
    )
  )

fig09 <- ggplot(df09, aes(x = raw_x, y = shap_x)) +
  annotate("rect", xmin = 200, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#8B4513", alpha = 0.06) +
  geom_vline(xintercept = 200, linetype = "dotted",
             color = "brown", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = elev_group, shape = elev_group),
             alpha = 0.55, size = 1.8) +
  geom_smooth(aes(color = elev_group, fill = elev_group),
              method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
  scale_color_manual(name = "Elevation\ngroup",
                     values = SHAP_GROUP_COLORS$elevation) +
  scale_fill_manual(name = "Elevation\ngroup",
                    values = SHAP_GROUP_COLORS$elevation) +
  scale_shape_manual(name = "Elevation\ngroup",
                     values = c("< 600 m a.s.l." = 16,
                                "\u2265 600 m a.s.l." = 17)) +
  annotate("text",
           x = 200, y = max(df09$shap_x, na.rm = TRUE) * 0.9,
           label = "Threshold\n(200 m)",
           size = 2.8, color = "brown",
           hjust = -0.1, fontface = "italic") +
  labs(
    x = "Sub-grid topographic complexity (SD of elevation, m)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig09_margin <- ggExtra::ggMarginal(fig09, type = "histogram", margins = "x",
                                     size = 8, fill = "gray80", color = "gray50")

save_figure(fig09_margin, fig_path("Fig09_SHAP_terrain.tiff"),
            width = 200, height = 150)

# =============================================================================
# Figure 10: Broadleaf fraction × phenology (leaf-on/leaf-off)
# =============================================================================
cat("[6/8] Figure 10: Broadleaf fraction x season...\n")

df10 <- extract_shap_pair(shap_long, frtp2_col, season_sum) %>%
  mutate(
    season_label = factor(
      ifelse(raw_color == 1,
             "Leaf-on (Summer)",
             "Leaf-off (Winter/Spring)"),
      levels = c("Leaf-off (Winter/Spring)", "Leaf-on (Summer)")
    )
  )

fig10 <- ggplot(df10, aes(x = raw_x, y = shap_x)) +
  geom_vline(xintercept = 0.3, linetype = "dotted",
             color = "darkgreen", linewidth = 0.6) +
  annotate("rect", xmin = 0.3, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#FFA500", alpha = 0.06) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = season_label, shape = season_label),
             alpha = 0.55, size = 1.8) +
  geom_smooth(aes(color = season_label, fill = season_label),
              method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
  scale_color_manual(name = "Phenological\nperiod",
                     values = SHAP_GROUP_COLORS$phenology) +
  scale_fill_manual(name = "Phenological\nperiod",
                    values = SHAP_GROUP_COLORS$phenology) +
  scale_shape_manual(name = "Phenological\nperiod",
                     values = c("Leaf-off (Winter/Spring)" = 17,
                                "Leaf-on (Summer)" = 16)) +
  annotate("text",
           x = 0.3, y = max(df10$shap_x, na.rm = TRUE) * 0.9,
           label = "Threshold\n(30%)",
           size = 2.8, color = "darkgreen",
           hjust = -0.1, fontface = "italic") +
  labs(
    x = "Broadleaf forest fraction (area ratio within footprint)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig10_margin <- ggExtra::ggMarginal(fig10, type = "histogram", margins = "x",
                                     size = 8, fill = "gray80", color = "gray50")

save_figure(fig10_margin, fig_path("Fig10_SHAP_phenology.tiff"),
            width = 200, height = 150)

# =============================================================================
# Figure 11: Distance to coast × eastness (foehn amplification zone)
# =============================================================================
cat("[7/8] Figure 11: Distance to coast x eastness...\n")

df11 <- extract_shap_pair(shap_long, dist_col, east_col) %>%
  mutate(
    dist_km = raw_x / 1000,
    aspect_group = factor(
      ifelse(raw_color > 0.3,
             "Leeward (East-facing)",
             "Other aspects"),
      levels = c("Other aspects", "Leeward (East-facing)")
    )
  )

fig11 <- ggplot(df11, aes(x = dist_km, y = shap_x)) +
  annotate("rect", xmin = 0, xmax = 25, ymin = -Inf, ymax = Inf,
           fill = "#FF6B6B", alpha = 0.08) +
  geom_vline(xintercept = 25, linetype = "dotted",
             color = "red3", linewidth = 0.6) +
  geom_vline(xintercept = 50, linetype = "dotted",
             color = "gray40", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = aspect_group), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(color = aspect_group, fill = aspect_group),
              method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.9) +
  scale_color_manual(name = "Aspect group",
                     values = SHAP_GROUP_COLORS$aspect) +
  scale_fill_manual(name = "Aspect group",
                    values = SHAP_GROUP_COLORS$aspect) +
  annotate("text", x = 12.5, y = max(df11$shap_x, na.rm = TRUE) * 0.85,
           label = "Foehn\namplification\nzone",
           size = 3, color = "red3",
           fontface = "italic", hjust = 0.5) +
  annotate("text", x = 25, y = min(df11$shap_x, na.rm = TRUE) * 0.85,
           label = "25 km", size = 2.8, color = "red3", hjust = -0.1) +
  annotate("text", x = 50, y = min(df11$shap_x, na.rm = TRUE) * 0.85,
           label = "50 km", size = 2.8, color = "gray40", hjust = -0.1) +
  labs(
    x = "Distance to coast (km)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig11_margin <- ggExtra::ggMarginal(fig11, type = "histogram", margins = "x",
                                     size = 8, fill = "gray80", color = "gray50")

save_figure(fig11_margin, fig_path("Fig11_SHAP_foehn.tiff"),
            width = 200, height = 150)

# =============================================================================
# Summary
# =============================================================================
cat("\n[8/8] Variable importance ranking (top 10):\n")

imp <- shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = round(mean(abs(value)), 4), .groups = "drop") %>%
  arrange(desc(mean_abs_shap))

print(head(imp, 10))

cat("\n================================================================\n")
cat("  07 Complete (manuscript figure numbering applied)\n")
cat("================================================================\n\n")
cat(sprintf("  Figure 6:  %s (SHAP summary)\n",
            fig_path("Fig06_SHAP_summary.tiff")))
cat(sprintf("  Figure 7:  %s (along-wind gradient)\n",
            fig_path("Fig07_SHAP_gradient.tiff")))
cat(sprintf("  Figure 8:  %s (downwind-upwind diff, labels fixed)\n",
            fig_path("Fig08_SHAP_updown.tiff")))
cat(sprintf("  Figure 9:  %s (terrain complexity)\n",
            fig_path("Fig09_SHAP_terrain.tiff")))
cat(sprintf("  Figure 10: %s (phenology)\n",
            fig_path("Fig10_SHAP_phenology.tiff")))
cat(sprintf("  Figure 11: %s (foehn)\n",
            fig_path("Fig11_SHAP_foehn.tiff")))

##############################################################################
# DIAGNOSTIC CODE — Run once before trusting Figure 8 orientation
##############################################################################
# 
# shap_long <- readRDS(inter_path("shap_long.rds"))
# tmp <- shap_long %>% filter(variable == "dem_updown_diff")
# 
# cat("rfvalue range:", range(tmp$rfvalue, na.rm = TRUE), "\n")
# cat("Mean SHAP for negative rfvalue:",
#     mean(tmp$value[tmp$rfvalue < 0], na.rm = TRUE), "\n")
# cat("Mean SHAP for positive rfvalue:",
#     mean(tmp$value[tmp$rfvalue > 0], na.rm = TRUE), "\n")
# 
# # Interpretation:
# #   Mean SHAP(negative) > Mean SHAP(positive) → rfvalue is "downwind - upwind"
# #     → Matches §2.6 definition. Keep Fig 8 code as-is.
# #
# #   Mean SHAP(negative) < Mean SHAP(positive) → rfvalue is "upwind - downwind"
# #     → Uncomment the line "raw_x = -raw_x" inside the Figure 8 block above.
##############################################################################
# =============================================================================
# Editable PPTX export (Fig 6–11)
#
# PowerPoint에서 범례 위치, 폰트 크기, 색상 등을 마우스로 직접 편집.
# 단, ggExtra::ggMarginal의 marginal histogram은 PPTX 편집 대상이 아니므로
# base plot만 export됨 (marginal histogram이 필요하면 기존 TIFF 사용).
# =============================================================================
cat("\n[+] Exporting editable PPTX...\n")

for (pkg in c("officer", "rvg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
}
library(officer)
library(rvg)

# 만들어진 figure만 모음 (컬럼 없어서 SKIP된 경우 자동 제외)
fig_list <- list()
if (exists("fig06")) fig_list[["Fig 06 — SHAP summary"]]         <- fig06
if (exists("fig07")) fig_list[["Fig 07 — Along-wind gradient"]]  <- fig07
if (exists("fig08")) fig_list[["Fig 08 — Downwind-upwind diff"]] <- fig08
if (exists("fig09")) fig_list[["Fig 09 — Terrain complexity"]]   <- fig09
if (exists("fig10")) fig_list[["Fig 10 — Broadleaf phenology"]]  <- fig10
if (exists("fig11")) fig_list[["Fig 11 — Coastal foehn"]]        <- fig11

# 16:9 슬라이드 한 장에 figure 하나씩
ppt <- read_pptx()

for (title in names(fig_list)) {
  ppt <- ppt %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = title,
            location = ph_location_type(type = "title")) %>%
    ph_with(value = dml(ggobj = fig_list[[title]]),
            location = ph_location(left = 0.5, top = 1.3,
                                   width = 12, height = 5.8))
}

pptx_path <- fig_path("SHAP_figures_editable.pptx")
print(ppt, target = pptx_path)

cat(sprintf("  → %s\n", pptx_path))
cat("  PowerPoint에서 편집 가능한 항목:\n")
cat("    · 범례 박스: 클릭 후 드래그하여 빈 공간으로 이동\n")
cat("    · 텍스트 (축 제목/눈금/범례 항목): 폰트 크기·색상 개별 변경\n")
cat("    · 점·선·영역: 그래픽 요소 개별 선택 후 속성 변경\n")