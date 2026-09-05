# ============================================================================
# 05_figures.R
# Publication figures for the two-arm study (reproduces the original style).
#
#   Figure 1 - Yearly NAION reporting rate per 10,000 drug reports for
#              semaglutide vs comparators (notoriety-bias fingerprint).
#   Figure 2 - Panel A: NAION disproportionality (ROR, log scale) by drug;
#              Panel B: drug-target MR (IVW) across ocular and control traits.
#
# Input : results/2_yearly_group_pt.csv    (from 01_faers_disproportionality.R)
#         results/2_yearly_group_total.csv (from 01_faers_disproportionality.R)
#         results/1_signal_table.csv       (from 01_faers_disproportionality.R)
#         results/B1_mr_results.csv        (from 04_drug_target_mr.R)
#         results/B1_leaveoneout.csv       (from 04_drug_target_mr.R)
# Output: figures/Figure1_NAION_yearly_trend.pdf / .png
#         figures/Figure2_two_arm_forest.pdf / .png
#         figures/FigureS1_leaveoneout.pdf / .png   (supplement)
#
# If the results CSVs are absent, the manuscript's key values (embedded below)
# are used so the figures can still be reproduced standalone.
# ============================================================================

library(data.table)
library(ggplot2)

OUT_DIR <- "results"
FIG_DIR <- "figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# NAION preferred terms (must match script 01)
pt_naion <- toupper(c(
  "Optic ischaemic neuropathy",
  "Non-arteritic anterior ischaemic optic neuropathy",
  "Anterior ischaemic optic neuropathy",
  "Ischaemic optic neuropathy"
))
pt_primary <- "OPTIC ISCHAEMIC NEUROPATHY"   # PT used for the primary ROR

# Drug-group order and colours (as in the published figures)
groups_plot <- c("Semaglutide", "Tirzepatide", "Liraglutide",
                 "SGLT2i", "DPP4i", "Metformin")
group_cols <- c(
  Semaglutide = "#C0392B", Tirzepatide = "#E67E22", Liraglutide = "#8064A2",
  SGLT2i = "#2E6F8E", DPP4i = "#5B9A5B", Metformin = "#7F8C8D"
)

theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey80"),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "plain", size = 14, hjust = 0.5)
  )

# ============================================================================
# Figure 1: yearly NAION reporting rate, semaglutide vs comparators
# ============================================================================

f_pt    <- file.path(OUT_DIR, "2_yearly_group_pt.csv")
f_total <- file.path(OUT_DIR, "2_yearly_group_total.csv")

if (file.exists(f_pt) && file.exists(f_total)) {
  pt_y  <- fread(f_pt);   setnames(pt_y,  tolower(names(pt_y)))
  tot_y <- fread(f_total); setnames(tot_y, tolower(names(tot_y)))

  cases <- pt_y[pt %chin% pt_naion, .(naion_cases = sum(a)), by = .(year, group)]
  yearly <- merge(cases,
                  tot_y[, .(year, group, reports = n_group)],
                  by = c("year", "group"), all = TRUE)
  yearly[is.na(naion_cases), naion_cases := 0]
  yearly[, rate := naion_cases / reports * 1e4]
  yearly <- yearly[group %chin% groups_plot]
} else {
  # Fallback: ILLUSTRATIVE values approximating the published figure.
  # For the final figure, always use the exported 2_yearly_*.csv files.
  message("2_yearly_*.csv not found - using embedded illustrative values.")
  yearly <- data.table(
    year = rep(2018:2026, times = 6),
    group = rep(groups_plot, each = 9),
    rate = c(
      0, 8, 3, 0, 3, 1, 51, 178, 202,      # Semaglutide
      rep(0, 6), 1, 12, 13,                # Tirzepatide
      3, 2, 0, 8, 1, 0, 9, 83, 14,         # Liraglutide
      rep(0, 8), 3,                        # SGLT2i
      0, 3, rep(0, 7),                     # DPP4i
      0, 0, 0, 0, 2, 15, 0, 0, 0           # Metformin
    )
  )
}
yearly[, group := factor(group, levels = groups_plot)]

p1 <- ggplot(yearly, aes(x = year, y = rate, colour = group)) +
  geom_line(data = ~ .x[group != "Semaglutide"], linewidth = 1.0) +
  geom_point(data = ~ .x[group != "Semaglutide"], size = 2.4) +
  geom_line(data = ~ .x[group == "Semaglutide"], linewidth = 1.8) +
  geom_point(data = ~ .x[group == "Semaglutide"], size = 3.2) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", colour = "black") +
  annotate("text", x = 2025.4, y = max(yearly$rate) * 0.14,
           label = "2026: partial year\n(Q1-Q2 only)",
           colour = "grey50", size = 3.6, hjust = 0.5) +
  scale_colour_manual(values = group_cols, breaks = groups_plot) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(linewidth = 2))) +
  scale_x_continuous(breaks = 2018:2026) +
  labs(x = "Year",
       y = "NAION reports per 10,000 drug reports") +
  theme_pub +
  theme(legend.key.width = unit(1.6, "lines"))

ggsave(file.path(FIG_DIR, "Figure1_NAION_yearly_trend.pdf"), p1,
       width = 9, height = 5.6)
ggsave(file.path(FIG_DIR, "Figure1_NAION_yearly_trend.png"), p1,
       width = 9, height = 5.6, dpi = 300)

# ============================================================================
# Figure 2, Panel A: NAION ROR by drug (primary PT, log scale)
# ============================================================================

f_signal <- file.path(OUT_DIR, "1_signal_table.csv")

if (file.exists(f_signal)) {
  sig <- fread(f_signal); setnames(sig, tolower(names(sig)))
  arm1 <- sig[pt == pt_primary & group != "GLP1RA_all",
              .(label = group, n = a, est = ror, lci = ror_lo, uci = ror_hi)]
  # keep the same drugs as the published panel
  keep_drugs <- c("Semaglutide", "Tirzepatide", "Liraglutide", "Dulaglutide",
                  "Exenatide", "Metformin", "DPP4i", "SGLT2i", "Atorvastatin")
  arm1 <- arm1[label %chin% keep_drugs]
  arm1[, label := factor(label, levels = rev(keep_drugs))]
} else {
  # Fallback: estimates shown in the published figure (primary PT)
  message("1_signal_table.csv not found - using embedded manuscript values.")
  arm1 <- data.table(
    label = factor(
      rev(c("Semaglutide", "Tirzepatide", "Liraglutide", "Dulaglutide",
            "Exenatide", "Metformin", "DPP4i", "SGLT2i", "Atorvastatin")),
      levels = rev(c("Semaglutide", "Tirzepatide", "Liraglutide", "Dulaglutide",
                     "Exenatide", "Metformin", "DPP4i", "SGLT2i", "Atorvastatin"))),
    n   = c(836, 134, 18, 22, 1, 11, 1, 4, 0),
    est = c(116.8, 7.5, 7.0, 4.5, 0.5, 1.6, 0.4, 0.3, 0.08),
    lci = c(106.8, 6.4, 4.4, 3.0, 0.07, 0.9, 0.05, 0.1, 0.01),
    uci = c(127.7, 8.8, 11.1, 6.8, 3.6, 2.9, 2.5, 1.0, 0.6)
  )
}
arm1[, ylab := paste0(label, "  (n=", n, ")")]
arm1[, ylab := factor(ylab, levels = ylab[match(levels(label), label)])]

p2a <- ggplot(arm1, aes(y = ylab, x = est)) +
  geom_vline(xintercept = 1, colour = "black", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0,
                 colour = "#C0392B", linewidth = 0.8) +
  geom_point(size = 3, colour = "#C0392B") +
  scale_x_log10(breaks = 10^(-2:3),
                labels = scales::math_format(10^.x),
                limits = c(0.008, 3000),
                minor_breaks = rep(2:9, 6) * 10^rep(-3:2, each = 8)) +
  annotation_logticks(sides = "b", outside = FALSE, linewidth = 0.3) +
  coord_cartesian(clip = "off") +
  labs(title = "A  Pharmacovigilance: NAION signal by drug",
       x = "ROR (95% CI), log scale - FAERS", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(size = 11),
        plot.margin = margin(5.5, 15, 5.5, 5.5))

# ============================================================================
# Figure 2, Panel B: drug-target MR (IVW) across outcomes
# ============================================================================

f_mr <- file.path(OUT_DIR, "B1_mr_results.csv")

mr_levels <- c("HbA1c (positive control)", "Intraocular pressure",
               "Vertical cup-disc ratio", "RNFL thickness", "GCIPL thickness",
               "Glaucoma (FinnGen)", "Skin tanning (negative control)")

oc_labels <- c(
  "finn-b-H7_GLAUCOMA"  = "Glaucoma (FinnGen)",
  "ebi-a-GCST004074"    = "Intraocular pressure",
  "ebi-a-GCST004075"    = "Vertical cup-disc ratio",
  "ebi-a-GCST90014266"  = "RNFL thickness",
  "ebi-a-GCST90014267"  = "GCIPL thickness",
  "ebi-a-GCST90014006"  = "HbA1c (positive control)",
  "ukb-b-533"           = "Skin tanning (negative control)"
)

if (file.exists(f_mr)) {
  mr_res <- fread(f_mr)
  arm2 <- mr_res[method == "Inverse variance weighted",
                 .(label = oc_labels[id.outcome], est = b,
                   lci = b - 1.96 * se, uci = b + 1.96 * se, p = pval)]
  arm2 <- arm2[!is.na(label)]
} else {
  # Fallback: IVW estimates reported in the manuscript
  message("B1_mr_results.csv not found - using embedded manuscript values.")
  arm2 <- data.table(
    label = c("Glaucoma (FinnGen)", "Intraocular pressure",
              "Vertical cup-disc ratio", "RNFL thickness", "GCIPL thickness",
              "HbA1c (positive control)", "Skin tanning (negative control)"),
    est = c( 0.167, -0.408, -0.014,  0.475,  1.060, -0.178,  0.024),
    se  = c( 0.241,  0.381,  0.018,  0.507,  0.667,  0.040,  0.028),
    p   = c(0.488, 0.284, 0.444, 0.349, 0.112, 7.6e-6, 0.402)
  )
  arm2[, `:=`(lci = est - 1.96 * se, uci = est + 1.96 * se)]
}
arm2[, label := factor(label, levels = rev(mr_levels))]
arm2[, kind := ifelse(grepl("control", label), "control", "ocular")]
arm2[, plab := ifelse(p < 0.001, sprintf("p = %.1e", p),
                      sprintf("p = %.2f", p))]

b_cols <- c(control = "#5B9A5B", ocular = "#2E6F8E")

p2b <- ggplot(arm2, aes(y = label, x = est)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lci, xmax = uci, colour = kind), height = 0,
                 linewidth = 0.8, show.legend = FALSE) +
  geom_point(aes(colour = kind), shape = 15, size = 3, show.legend = FALSE) +
  geom_text(aes(x = 3.15, label = plab), colour = "grey55", size = 3.6,
            hjust = 1) +
  scale_colour_manual(values = b_cols) +
  scale_x_continuous(breaks = -2:3, limits = c(-2.2, 3.2)) +
  coord_cartesian(clip = "off") +
  labs(title = "B  Drug-target MR across control and ocular outcomes",
       x = "IVW estimate per 1-SD GLP1R expression (95% CI)", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(size = 11),
        plot.margin = margin(5.5, 15, 5.5, 5.5))

# ---- combine and save --------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p2 <- p2a + p2b + plot_layout(widths = c(1, 1.05))
  ggsave(file.path(FIG_DIR, "Figure2_two_arm_forest.pdf"), p2,
         width = 13, height = 6)
  ggsave(file.path(FIG_DIR, "Figure2_two_arm_forest.png"), p2,
         width = 13, height = 6, dpi = 300)
} else {
  message("patchwork not installed - saving panels A and B separately.")
  ggsave(file.path(FIG_DIR, "Figure2A_faers_forest.pdf"), p2a, width = 6.5, height = 6)
  ggsave(file.path(FIG_DIR, "Figure2A_faers_forest.png"), p2a, width = 6.5, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, "Figure2B_mr_forest.pdf"), p2b, width = 7, height = 6)
  ggsave(file.path(FIG_DIR, "Figure2B_mr_forest.png"), p2b, width = 7, height = 6, dpi = 300)
}

# ============================================================================
# Figure S1 (supplement): leave-one-out analysis across ocular outcomes
# ============================================================================
# Input: results/B1_leaveoneout.csv (from 04_drug_target_mr.R, mr_leaveoneout).
# Each row re-estimates the IVW effect after removing one instrument; the
# "All" row is the overall IVW estimate (red square).

f_loo <- file.path(OUT_DIR, "B1_leaveoneout.csv")

if (file.exists(f_loo)) {
  loo <- fread(f_loo)
  loo[, oc := oc_labels[id.outcome]]
  loo <- loo[!is.na(oc)]

  loo_plot <- unique(loo[, .(oc, SNP, b, se)], by = c("oc", "SNP"))
  loo_plot[, is_all := grepl("All", SNP)]
  loo_plot[, snp_lab := ifelse(is_all, "All (IVW)", SNP)]
  # one common y order for all panels (same 6 instruments everywhere):
  # "All (IVW)" at the bottom, SNPs above it sorted by mean beta
  snp_ord <- loo_plot[!is_all == TRUE, .(mb = mean(b)), by = snp_lab][order(mb)]$snp_lab
  loo_plot[, snp_lab := factor(snp_lab, levels = c("All (IVW)", snp_ord))]

  pS1 <- ggplot(loo_plot, aes(x = b, y = snp_lab)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se,
                       colour = is_all),
                   height = 0, linewidth = 0.8, show.legend = FALSE) +
    geom_point(aes(colour = is_all, shape = is_all), size = 2.4,
               show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "#2E6F8E", `TRUE` = "#C0392B")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
    facet_wrap(~ oc, scales = "free_x", ncol = 3) +
    labs(title = "Leave-one-out analysis: GLP1R drug-target MR",
         x = "IVW estimate leaving out each SNP (95% CI)", y = NULL) +
    theme_pub +
    theme(strip.text = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 8))

  ggsave(file.path(FIG_DIR, "FigureS1_leaveoneout.pdf"), pS1,
         width = 9.5, height = 6.5)
  ggsave(file.path(FIG_DIR, "FigureS1_leaveoneout.png"), pS1,
         width = 9.5, height = 6.5, dpi = 300)
} else {
  message("B1_leaveoneout.csv not found - skipping Figure S1.")
}
