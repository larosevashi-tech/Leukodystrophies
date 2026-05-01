# ============================================================
# BSA Standard Curve Analysis
# ============================================================

# ── 1. Libraries ─────────────────────────────────────────────
library(readxl)
library(ggplot2)
library(dplyr)
library(broom)
library(writexl)
library(openxlsx)

# ── 2. Load & prep data ───────────────────────────────────────
data <- read_excel(r"(C:\Users\Bernard420\OneDrive\Documents\Lab_Projects\2026-04-24_BSA_Analysis\BSA_Data.xlsx)")

# Mean_Abs — average absorbance across replicates, ignoring NAs
# SD_Abs — standard deviation of absorbance across replicates
# CV_pct — coefficient of variation (SD/Mean × 100), expressing variability as a percentage 
#         — useful for checking reproducibility between replicates
# n_reps — count of how many replicate rows were in that group

data_avg <- data %>%
  group_by(sample_id, conc_mg_mL) %>%
  summarise(
    Mean_Abs = mean(absorbance, na.rm = TRUE),
    SD_Abs   = sd(absorbance,   na.rm = TRUE),
    CV_pct   = (SD_Abs / Mean_Abs) * 100,
    n_reps   = n(),
    .groups  = "drop"
  )

# blank correction
blank <- data_avg %>% filter(sample_id == "blank") %>% pull(Mean_Abs)
data_avg <- data_avg %>% mutate(Mean_Abs = Mean_Abs - blank)

standards <- data_avg %>%
  filter(!is.na(conc_mg_mL), CV_pct <= 20 | n_reps == 1) %>%
  filter(!sample_id %in% c("blank", "empty"))

unknowns <- data_avg %>% filter(is.na(conc_mg_mL))


# ── 3. Fit Linear Model (Standard Curve) ─────────────────────
model         <- lm(Mean_Abs ~ conc_mg_mL, data = standards)
model_summary <- summary(model)
intercept     <- coef(model)[1]
slope         <- coef(model)[2]
r_squared_lin <- model_summary$r.squared

cat("─── BSA Standard Curve: Linear ─────────────────────────\n")
cat(sprintf("  Mean_Abs = %.5f + %.5f × [Protein] (mg/mL)\n", intercept, slope))
cat(sprintf("  R²       = %.4f\n", r_squared_lin))
cat("────────────────────────────────────────────────────────\n\n")

# ── 4. Fit Quadratic Model (Standard Curve) ───────────────────
model_poly    <- lm(Mean_Abs ~ poly(conc_mg_mL, 2, raw = TRUE), data = standards)
coefs         <- coef(model_poly)
r_squared_poly <- summary(model_poly)$r.squared

cat("─── BSA Standard Curve: Quadratic ──────────────────────\n")
cat(sprintf("  Mean_Abs = %.5f + %.5fx + %.5fx²\n", coefs[1], coefs[2], coefs[3]))
cat(sprintf("  R²       = %.4f\n", r_squared_poly))
cat("────────────────────────────────────────────────────────\n\n")

# ── 5. Compare Models ─────────────────────────────────────────
cat("─── Model Comparison ────────────────────────────────────\n")
cat(sprintf("  Linear R²    = %.4f\n", r_squared_lin))
cat(sprintf("  Quadratic R² = %.4f\n", r_squared_poly))
if (r_squared_poly > r_squared_lin) {
  cat("  → Quadratic fit is better\n")
} else {
  cat("  → Linear fit is better\n")
}
cat("────────────────────────────────────────────────────────\n\n")

# ── 6. Smooth curve for plotting ──────────────────────────────
curve_range <- data.frame(
  conc_mg_mL = seq(min(standards$conc_mg_mL),
                   max(standards$conc_mg_mL), length.out = 300)
)
curve_range$fitted_abs <- coefs[1] +
  coefs[2] * curve_range$conc_mg_mL +
  coefs[3] * curve_range$conc_mg_mL^2

# ── 7. Standard recovery, Back-calculate unknowns ─────────────────────────
if (r_squared_poly > r_squared_lin) {
  
  standards <- standards %>%
    mutate(
      fitted_abs     = coefs[1] + coefs[2]*conc_mg_mL + coefs[3]*conc_mg_mL^2,
      residual       = Mean_Abs - fitted_abs,
      c_val          = coefs[1] - Mean_Abs,
      discriminant   = coefs[2]^2 - 4 * coefs[3] * c_val,
      conc_recovered = (-coefs[2] + sqrt(discriminant)) / (2 * coefs[3]),
      recovery_pct   = if_else(conc_mg_mL > 0,
                               (conc_recovered / conc_mg_mL) * 100, NA_real_)
    ) %>%
    select(-c_val, -discriminant)
  
  unknowns <- unknowns %>%
    mutate(
      c_val         = coefs[1] - Mean_Abs,
      discriminant  = coefs[2]^2 - 4 * coefs[3] * c_val,
      conc_mg_mL    = (-coefs[2] + sqrt(discriminant)) / (2 * coefs[3]),
      abs_predicted = coefs[1] + coefs[2]*conc_mg_mL + coefs[3]*conc_mg_mL^2,
      abs_residual  = Mean_Abs - abs_predicted,
      in_range      = case_when(
        conc_mg_mL < min(standards$conc_mg_mL) ~ "Below range",
        conc_mg_mL > max(standards$conc_mg_mL) ~ "Above range",
        TRUE                                    ~ "In range"
      )
    ) %>%
    select(-c_val, -discriminant)

} else {
  
  standards <- standards %>%
    mutate(
      fitted_abs     = intercept + slope * conc_mg_mL,
      residual       = Mean_Abs - fitted_abs,
      conc_recovered = (Mean_Abs - intercept) / slope,
      recovery_pct   = if_else(conc_mg_mL > 0,
                               (conc_recovered / conc_mg_mL) * 100, NA_real_)
    )
  
  unknowns <- unknowns %>%
    mutate(
      conc_mg_mL    = (Mean_Abs - intercept) / slope,
      abs_predicted = intercept + slope * conc_mg_mL,
      abs_residual  = Mean_Abs - abs_predicted,
      in_range      = case_when(
        conc_mg_mL < min(standards$conc_mg_mL) ~ "Below range",
        conc_mg_mL > max(standards$conc_mg_mL) ~ "Above range",
        TRUE                                    ~ "In range"
      )
    )
}
# ════════════════════════════════════════════════════════════
# AUTO-SELECT MODEL
# ════════════════════════════════════════════════════════════
use_quadratic <- r_squared_poly > r_squared_lin
cat(sprintf("Using %s model for plots\n", ifelse(use_quadratic, "quadratic", "linear")))

if (use_quadratic) {
  curve_range$fitted_abs <- coefs[1] +
    coefs[2] * curve_range$conc_mg_mL +
    coefs[3] * curve_range$conc_mg_mL^2
  eq_label <- sprintf("y = %.4fx² + %.4fx + %.4f\nR² = %.4f",
                      coefs[3], coefs[2], coefs[1], r_squared_poly)
  plot_subtitle_model <- "Quadratic Fit"
} else {
  curve_range$fitted_abs <- intercept + slope * curve_range$conc_mg_mL
  eq_label <- sprintf("y = %.4f + %.4fx\nR² = %.4f",
                      intercept, slope, r_squared_lin)
  plot_subtitle_model <- "Linear Fit"
}

# ════════════════════════════════════════════════════════════
# PLOTS
# ════════════════════════════════════════════════════════════

# ── Shared theme ──────────────────────────────────────────────
bsa_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", size = 15, color = "#1e3a5f"),
      plot.subtitle   = element_text(color = "grey40", size = 10),
      axis.title      = element_text(color = "#1e3a5f"),
      axis.text       = element_text(color = "grey30"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4)
    )
}

BLUE  <- "#1D4ED8"
RED   <- "#EF4444"
DARK  <- "#1e3a5f"


# ── Plot 1: Standards curve with error bars ───────────────────
p_curve <- ggplot() +
  geom_ribbon(data = curve_range,
              aes(x = conc_mg_mL,
                  ymin = fitted_abs - 0.01,
                  ymax = fitted_abs + 0.01),
              fill = BLUE, alpha = 0.10) +
  geom_line(data = curve_range,
            aes(x = conc_mg_mL, y = fitted_abs),
            color = BLUE, linewidth = 1) +
  geom_errorbar(data = standards,
                aes(x = conc_mg_mL,
                    ymin = Mean_Abs - SD_Abs,
                    ymax = Mean_Abs + SD_Abs),
                width = 0.03, color = BLUE, alpha = 0.6) +
  geom_point(data = standards,
             aes(x = conc_mg_mL, y = Mean_Abs),
             color = BLUE, size = 1, shape = 16) +
  annotate("text",
           x = min(standards$conc_mg_mL),
           y = max(standards$Mean_Abs) * 0.97,
           label = eq_label,
           hjust = 0, size = 3.5, fontface = "italic", color = DARK) +
  labs(title    = paste("BSA Standard Curve —", plot_subtitle_model),
       subtitle = "Error bars = SD across replicates",
       x        = "Protein Concentration (mg/mL)",
       y        = "Absorbance (562 nm)") +
  bsa_theme()

# ── Plot 2: Unknowns projected onto curve with error bars ─────
p_unknowns <- ggplot() +
  geom_ribbon(data = curve_range,
              aes(x = conc_mg_mL,
                  ymin = fitted_abs - 0.01,
                  ymax = fitted_abs + 0.01),
              fill = BLUE, alpha = 0.10) +
  geom_line(data = curve_range,
            aes(x = conc_mg_mL, y = fitted_abs),
            color = BLUE, linewidth = 1) +
  geom_errorbar(data = standards,
                aes(x = conc_mg_mL,
                    ymin = Mean_Abs - SD_Abs,
                    ymax = Mean_Abs + SD_Abs),
                width = 0.03, color = BLUE, alpha = 0.4) +
  geom_point(data = standards,
             aes(x = conc_mg_mL, y = Mean_Abs),
             color = BLUE, size = 1, shape = 16) +
  geom_segment(data = unknowns,
               aes(x = conc_mg_mL, xend = conc_mg_mL,
                   y = -Inf, yend = Mean_Abs),
               linetype = "dotted", color = RED, alpha = 0.6) +
  geom_segment(data = unknowns,
               aes(x = -Inf, xend = conc_mg_mL,
                   y = Mean_Abs, yend = Mean_Abs),
               linetype = "dotted", color = RED, alpha = 0.6) +
  geom_errorbar(data = unknowns,
                aes(x = conc_mg_mL,
                    ymin = Mean_Abs - SD_Abs,
                    ymax = Mean_Abs + SD_Abs),
                width = 0.03, color = RED, alpha = 0.6) +
  geom_point(data = unknowns,
             aes(x = conc_mg_mL, y = Mean_Abs),
             color = RED, size = 1, shape = 17) +
  geom_text(data = unknowns,
            aes(x = conc_mg_mL, y = Mean_Abs, label = sample_id),
            vjust = -1, size = 3.2, color = RED) +
  annotate("text",
           x = min(standards$conc_mg_mL),
           y = max(standards$Mean_Abs) * 0.97,
           label = eq_label,
           hjust = 0, size = 3.5, fontface = "italic", color = DARK) +
  labs(title    = paste("BSA Standard Curve —", plot_subtitle_model, "| Unknowns Projected"),
       subtitle = "Blue circles = standards  |  Red triangles = unknowns  |  Error bars = SD",
       x        = "Protein Concentration (mg/mL)",
       y        = "Absorbance (562 nm)") +
  bsa_theme()

# ── Plot 3: Residual plot ─────────────────────────────────────
p_residual <- ggplot(standards, aes(x = conc_mg_mL, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_segment(aes(xend = conc_mg_mL, yend = 0),
               color = BLUE, alpha = 0.4, linewidth = 0.7) +
  geom_point(color = BLUE, size = 3.5, shape = 16) +
  geom_text(aes(label = sample_id),
            vjust = -1, size = 3, color = DARK) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  labs(title    = paste("Residual Plot —", plot_subtitle_model),
       subtitle = "Residual = Observed − Fitted  |  Points near 0 = good fit",
       x        = "Protein Concentration (mg/mL)",
       y        = "Residual (Absorbance units)") +
  bsa_theme()

# ── Plot 4: CV% per standard ──────────────────────────────────
p_cv <- ggplot(standards, aes(x = factor(conc_mg_mL), y = CV_pct)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = RED, alpha = 0.6) +
  geom_col(fill = BLUE, alpha = 0.7, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", CV_pct)),
            vjust = -0.5, size = 3.2, color = DARK) +
  annotate("text", x = Inf, y = 20,
           label = "CV threshold", hjust = 1.1, vjust = -0.5,
           size = 3, color = RED) +
  labs(title    = "CV% per Standard Concentration",
       subtitle = "Dashed line = 20% threshold",
       x        = "Protein Concentration (mg/mL)",
       y        = "CV (%)") +
  bsa_theme()

# ── Plot 5: Recovery % per standard ──────────────────────────
p_recovery <- ggplot(standards %>% filter(!is.na(recovery_pct)),
                     aes(x = factor(conc_mg_mL), y = recovery_pct)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_col(fill = BLUE, alpha = 0.7, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", recovery_pct)),
            vjust = -0.5, size = 3.2, color = DARK) +
  labs(title    = "Standard Recovery %",
       subtitle = "Ideal recovery = 100%  |  Back-calculated vs known concentration",
       x        = "Protein Concentration (mg/mL)",
       y        = "Recovery (%)") +
  bsa_theme()

print(p_curve)
print(p_unknowns)
print(p_residual)
print(p_cv)
print(p_recovery)

# ════════════════════════════════════════════════════════════
# EXPORT TO EXCEL
# ════════════════════════════════════════════════════════════

insert_plot <- function(wb, sheet, plot_obj, tmpfile,
                        row = 2, col = 2, w = 7, h = 5) {
  ggsave(tmpfile, plot = plot_obj, width = w, height = h, dpi = 300)
  insertImage(wb, sheet, tmpfile,
              startRow = row, startCol = col,
              width = w, height = h, units = "in")
  file.remove(tmpfile)
}

wb <- createWorkbook()

addWorksheet(wb, "Raw Data")
writeDataTable(wb, "Raw Data",
               data,
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Averaged Data")
writeDataTable(wb, "Averaged Data",
               data_avg %>% select(sample_id, conc_mg_mL,
                                   Mean_Abs, SD_Abs, CV_pct, n_reps),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Standards")
writeDataTable(wb, "Standards",
               standards %>% select(sample_id, conc_mg_mL, Mean_Abs, SD_Abs,
                                    fitted_abs, residual,
                                    conc_recovered, recovery_pct),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Unknowns")
writeDataTable(wb, "Unknowns",
               unknowns %>% select(sample_id, Mean_Abs, SD_Abs, conc_mg_mL,
                                   abs_predicted, abs_residual, in_range),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Standard Curve")
insert_plot(wb, "Standard Curve",   p_curve,    "tmp_curve.png")

addWorksheet(wb, "Curve + Unknowns")
insert_plot(wb, "Curve + Unknowns", p_unknowns, "tmp_unknowns.png")

addWorksheet(wb, "Residual Plot")
insert_plot(wb, "Residual Plot",    p_residual, "tmp_residual.png")

addWorksheet(wb, "CV Plot")
insert_plot(wb, "CV Plot",          p_cv,       "tmp_cv.png")

addWorksheet(wb, "Recovery Plot")
insert_plot(wb, "Recovery Plot",    p_recovery, "tmp_recovery.png")

saveWorkbook(
  wb,
  "C:/Users/Bernard420/OneDrive/Documents/Lab_Projects/2026-04-24_BSA_Analysis/BSA_Full_Results.xlsx",
  overwrite = TRUE
)
cat("Done → BSA_Full_Results.xlsx\n")
