## ============================================================================
## Arm 2 — Instrument selection
## GLP1R cis-eQTL instruments from eQTLGen (whole blood, FDR < 0.05)
##
## Input: eQTLGen significant cis-eQTL file, downloaded from
##        https://eqtlgen.org/cis-eqtls.html
##        (2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz)
## Output: glp1r_lead_snps.csv (LD-independent lead instruments)
##
## Note: ld_clump() queries the OpenGWAS API and requires a JWT token;
##       see https://api.opengwas.io/ (free registration via GitHub).
## ============================================================================
library(data.table)
library(ieugwasr)   # if needed: remotes::install_github("MRCIEU/ieugwasr")

## ---------- Configuration (EDIT THESE PATHS) ----------
EQTL_FILE <- "/path/to/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"  ##Download from https://eqtlgen.org/cis-eqtls.html
OUT_FILE  <- "/path/to/glp1r_lead_snps.csv"

## ---------- 1. Extract GLP1R cis-eQTLs ----------
glp1r <- fread(EQTL_FILE, header = TRUE)
glp1r <- glp1r[glp1r$GeneSymbol == "GLP1R", ]
cat("Significant cis-eQTLs for GLP1R:", nrow(glp1r), "\n")

## ---------- 2. LD clumping (1000 Genomes EUR, r2 < 0.1, 500 kb window) ----------
clump_input <- glp1r[, .(rsid = SNP, pval = as.numeric(Pvalue))]
clumped <- ld_clump(clump_input,
                    clump_r2 = 0.1,
                    clump_kb = 500,
                    pop = "EUR")
cat("Independent SNPs after clumping:", nrow(clumped), "\n")
print(clumped)

## ---------- 3. Export lead instruments ----------
lead <- glp1r[SNP %chin% clumped$rsid,
              .(SNP, SNPPos, AssessedAllele, OtherAllele, Zscore, Pvalue, NrSamples)]
lead[, absZ := abs(as.numeric(Zscore))]
setorder(lead, -absZ)
print(lead)
fwrite(lead, OUT_FILE)
