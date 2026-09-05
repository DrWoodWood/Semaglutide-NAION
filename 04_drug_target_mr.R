# ============================================================================
# 04_drug_target_mr.R
# Arm 2: Drug-target Mendelian randomization of GLP1R expression against
#        glaucoma-spectrum ocular outcomes (proxies for NAION), with positive
#        (HbA1c) and negative (skin tanning) controls.
#
# Input : results/glp1r_lead_snps.csv  (from 03_glp1r_instruments.R)
# Output: results/B1_mr_results.csv
#         results/B1_heterogeneity.csv
#         results/B1_pleiotropy.csv
#         results/B1_leaveoneout.csv
#         results/B1_mr_results_noLowN.csv   (sensitivity, excl. rs114977861)
#
# Requires an OpenGWAS JWT token: Sys.setenv(OPENGWAS_JWT = "<your token>")
# ============================================================================

library(data.table)
library(TwoSampleMR)
library(ieugwasr)

# ---- 0. Configuration -------------------------------------------------------
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LEAD_FILE <- file.path(OUT_DIR, "glp1r_lead_snps.csv")

# ---- 1. Exposure: GLP1R cis-eQTL lead SNPs ---------------------------------
lead <- fread(LEAD_FILE)

# Convert eQTLGen Z-scores to beta/se (approximation, N per SNP)
lead[, beta := as.numeric(Zscore) / sqrt(as.numeric(NrSamples))]
lead[, se   := 1 / sqrt(as.numeric(NrSamples))]

# Note: the eQTLGen significant cis-eQTL file has no MAF column, so no
# eaf_col is passed. Palindromic SNPs with unknown MAF will be flagged as
# ambiguous at harmonisation (action = 2); check the harmonise output and
# supply EUR MAF manually only if a lead SNP is dropped for that reason.
exp_dat <- format_data(
  as.data.frame(lead),
  type            = "exposure",
  snp_col         = "SNP",
  beta_col        = "beta",
  se_col          = "se",
  effect_allele_col  = "AssessedAllele",
  other_allele_col   = "OtherAllele",
  pval_col        = "Pvalue"
)
exp_dat$exposure <- "GLP1R expression (eQTLGen whole blood)"

# ---- 2. Outcomes -------------------------------------------------------------
# Glaucoma-spectrum traits used as NAION proxies:
#   disc at risk (vCDR), optic nerve head perfusion (IOP),
#   retinal ganglion cell integrity (RNFL, GCIPL), glaucoma liability.
outcome_ids <- c(
  glaucoma = "finn-b-H7_GLAUCOMA",
  IOP      = "ebi-a-GCST004074",
  vCDR     = "ebi-a-GCST004075",
  rnfl     = "ebi-a-GCST90014266",
  gcipl    = "ebi-a-GCST90014267",
  hba1c    = "ebi-a-GCST90014006",  # positive control
  tanning  = "ukb-b-533"            # negative control
)

# sanity check that all outcome datasets are accessible
info <- gwasinfo(outcome_ids)
print(as.data.frame(info)[, intersect(
  c("id", "trait", "sample_size", "nsnp", "year", "population"),
  names(info))])

# ---- 3. Extract, harmonise, run MR ------------------------------------------
out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_ids)
dat     <- harmonise_data(exp_dat, out_dat, action = 2)

res <- mr(
  dat,
  method_list = c("mr_ivw", "mr_egger_regression",
                  "mr_weighted_median", "mr_wald_ratio")
)
het <- mr_heterogeneity(dat)
ple <- mr_pleiotropy_test(dat)
loo <- mr_leaveoneout(dat)

fwrite(res, file.path(OUT_DIR, "B1_mr_results.csv"))
fwrite(het, file.path(OUT_DIR, "B1_heterogeneity.csv"))
fwrite(ple, file.path(OUT_DIR, "B1_pleiotropy.csv"))
fwrite(loo, file.path(OUT_DIR, "B1_leaveoneout.csv"))

# ---- 4. Sensitivity analysis -------------------------------------------------
# rs114977861 has a much smaller instrument sample size (NrSamples = 3,243);
# repeat the analysis without it.
dat_sens <- dat[dat$SNP != "rs114977861", ]

res_sens <- mr(
  dat_sens,
  method_list = c("mr_ivw", "mr_egger_regression",
                  "mr_weighted_median", "mr_wald_ratio")
)
fwrite(res_sens, file.path(OUT_DIR, "B1_mr_results_noLowN.csv"))
