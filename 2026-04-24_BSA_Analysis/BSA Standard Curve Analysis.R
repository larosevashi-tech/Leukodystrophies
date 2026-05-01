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

# ── 2. Load & prep data ───────────────────────────────────────
data <- read_excel(r"(C:\Users\Bernard420\OneDrive\Documents\Lab_Projects\2026-04-24_BSA_Analysis\BSA_Data.xlsx)")

data_avg <- data %>%
  group_by(sample_id, conc_mg_mL) %>%
  summarise(
    Mean_Abs = mean(absorbance, na.rm = TRUE),
    SD_Abs   = sd(absorbance,   na.rm = TRUE),
    CV_pct   = (SD_Abs / Mean_Abs) * 100,
    n_reps   = n(),
    .groups  = "drop"
  )

# RIPA blank correction
ripa_blank <- data_avg %>% filter(sample_id == "RIPA") %>% pull(Mean_Abs)
data_avg   <- data_avg %>% mutate(Mean_Abs = Mean_Abs - ripa_blank)

standards <- data_avg %>%
  filter(!is.na(conc_mg_mL), sample_id != "BSA_0.125")
unknowns  <- data_avg %>% filter(is.na(conc_mg_mL))

# ── 3. Fit quadratic model ────────────────────────────────────
model_poly <- lm(Mean_Abs ~ poly(conc_mg_mL, 2, raw = TRUE), data = standards)
coefs      <- coef(model_poly)
r_squared  <- summary(model_poly)$r.squared

eq_label <- sprintf("y = %.4fx² + %.4fx + %.4f\nR² = %.4f",
                    coefs[3], coefs[2], coefs[1], r_squared)

curve_range <- data.frame(
  conc_mg_mL = seq(min(standards$conc_mg_mL),
                   max(standards$conc_mg_mL), length.out = 300)
)
curve_range$fitted_abs <- coefs[1] +
  coefs[2] * curve_range$conc_mg_mL +
  coefs[3] * curve_range$conc_mg_mL^2

# ── 4. Standard recovery & residuals ─────────────────────────
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

# ── 5. Back-calculate unknowns ────────────────────────────────
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


# ════════════════════════════════════════════════════════════
# PLOTS
# ════════════════════════════════════════════════════════════

# ── Plot 1: Standards-only curve ─────────────────────────────
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
             color = BLUE, size = 3, shape = 16) +
  annotate("text",
           x = min(standards$conc_mg_mL),
           y = max(standards$Mean_Abs) * 0.97,
           label = eq_label,
           hjust = 0, size = 3.5, fontface = "italic", color = DARK) +
  labs(title    = "BSA Standard Curve — Quadratic Fit",
       subtitle = "Error bars = SD across replicates",
       x        = "Protein Concentration (mg/mL)",
       y        = "Absorbance (562 nm)") +
  bsa_theme()

# ── Plot 2: Curve + unknowns projected ───────────────────────
p_unknowns <- ggplot() +
  geom_ribbon(data = curve_range,
              aes(x = conc_mg_mL,
                  ymin = fitted_abs - 0.01,
                  ymax = fitted_abs + 0.01),
              fill = BLUE, alpha = 0.10) +
  geom_line(data = curve_range,
            aes(x = conc_mg_mL, y = fitted_abs),
            color = BLUE, linewidth = 1) +
  geom_point(data = standards,
             aes(x = conc_mg_mL, y = Mean_Abs),
             color = BLUE, size = 3, shape = 16) +
  # Crosshairs
  geom_segment(data = unknowns,
               aes(x = conc_mg_mL, xend = conc_mg_mL,
                   y = -Inf, yend = Mean_Abs),
               linetype = "dotted", color = RED, alpha = 0.6) +
  geom_segment(data = unknowns,
               aes(x = -Inf, xend = conc_mg_mL,
                   y = Mean_Abs, yend = Mean_Abs),
               linetype = "dotted", color = RED, alpha = 0.6) +
  geom_point(data = unknowns,
             aes(x = conc_mg_mL, y = Mean_Abs),
             color = RED, size = 4, shape = 17) +
  geom_text(data = unknowns,
            aes(x = conc_mg_mL, y = Mean_Abs, label = sample_id),
            vjust = -1, size = 3.2, color = RED) +
  annotate("text",
           x = min(standards$conc_mg_mL),
           y = max(standards$Mean_Abs) * 0.97,
           label = eq_label,
           hjust = 0, size = 3.5, fontface = "italic", color = DARK) +
  labs(title    = "BSA Standard Curve — Unknowns Projected",
       subtitle = "Blue circles = standards  |  Red triangles = unknowns",
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
  labs(title    = "Residual Plot — Quadratic Fit",
       subtitle = "Residual = Observed − Fitted  |  Points near 0 = good fit",
       x        = "Protein Concentration (mg/mL)",
       y        = "Residual (Absorbance units)") +
  bsa_theme()

print(p_curve)
print(p_unknowns)
print(p_residual)


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

addWorksheet(wb, "Averaged Data")
writeDataTable(wb, "Averaged Data",
               data_avg %>% select(sample_id, conc_mg_mL,
                                   Mean_Abs, SD_Abs, CV_pct, n_reps),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Standards")
writeDataTable(wb, "Standards",
               standards %>% select(sample_id, conc_mg_mL, Mean_Abs,
                                    fitted_abs, residual,
                                    conc_recovered, recovery_pct),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Unknowns")
writeDataTable(wb, "Unknowns",
               unknowns %>% select(sample_id, Mean_Abs, conc_mg_mL,
                                   abs_predicted, abs_residual, in_range),
               tableStyle = "TableStyleMedium9")

addWorksheet(wb, "Standard Curve")
insert_plot(wb, "Standard Curve",   p_curve,    "tmp_curve.png")

addWorksheet(wb, "Curve + Unknowns")
insert_plot(wb, "Curve + Unknowns", p_unknowns, "tmp_unknowns.png")

addWorksheet(wb, "Residual Plot")
insert_plot(wb, "Residual Plot",    p_residual, "tmp_residual.png")

saveWorkbook(
  wb,
  "C:/Users/Bernard420/OneDrive/Documents/Lab_Projects/2026-04-24_BSA_Analysis/BSA_Full_Results.xlsx",
  overwrite = TRUE
)
cat("Done → BSA_Full_Results.xlsx\n")