## ============================================================================
## Arm 1 — Supplementary analysis (lightweight, per-quarter streaming)
## Indication composition + brand stratification of NAION x semaglutide cases
##
## Memory-safe design: reads one quarterly slim RDS at a time, filters early,
## never holds the full drug table in memory.
##
## Usage: run AFTER 01_faers_disproportionality.R has built rds_cache.
##        Restart R first (RStudio: Session -> Restart R), then source this file.
## ============================================================================
library(data.table)

## ---------- Configuration (EDIT THESE PATHS) ----------
FAERS_DIR <- "/path/to/faers"
CACHE_DIR <- file.path(FAERS_DIR, "rds_cache")
OUT_DIR   <- file.path(FAERS_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(CACHE_DIR, "\\.rds$", full.names = TRUE)
if (!length(files)) stop("No cache files in rds_cache; run 01_faers_disproportionality.R first")

pt_naion <- c("OPTIC ISCHAEMIC NEUROPATHY",
              "NON-ARTERITIC ANTERIOR ISCHAEMIC OPTIC NEUROPATHY",
              "ANTERIOR ISCHAEMIC OPTIC NEUROPATHY",
              "ISCHAEMIC OPTIC NEUROPATHY")
sema_rx <- "SEMAGLUTIDE|OZEMPIC|WEGOVY|RYBELSUS"

res_indi  <- list()
res_brand <- list()

for (f in files) {
  tag <- basename(f)
  x <- readRDS(f)
  demo <- x$demo; drug <- x$drug; reac <- x$reac; indi <- x$indi
  rm(x); gc()

  ## Within-quarter dedup (latest version only)
  demo[, caseversion := as.integer(caseversion)]
  setorder(demo, caseid, -caseversion)
  keep <- demo[!duplicated(caseid), primaryid]

  ## Filter NAION cases first (very few), then semaglutide cases —
  ## avoids materializing database-wide string columns
  naion_ids <- unique(reac[pt %chin% pt_naion & primaryid %chin% keep, primaryid])
  if (!length(naion_ids)) next
  sema_ids <- drug[role_cod == "PS" &
                     (grepl(sema_rx, toupper(drugname)) |
                      grepl("SEMAGLUTIDE", toupper(prod_ai))), primaryid]
  hit <- intersect(naion_ids, intersect(sema_ids, keep))
  if (!length(hit)) next

  ## Indication records
  res_indi[[tag]] <- indi[primaryid %chin% hit, .(primaryid, indi_pt)]

  ## Brand stratification
  res_brand[[tag]] <- unique(drug[primaryid %chin% hit & role_cod == "PS",
                                  .(primaryid, dn = toupper(drugname))])[
    , .(Ozempic  = any(grepl("OZEMPIC",  dn)),
        Wegovy   = any(grepl("WEGOVY",   dn)),
        Rybelsus = any(grepl("RYBELSUS", dn))),
    by = primaryid]

  rm(demo, drug, reac, indi); gc()
  cat(tag, " NAION x semaglutide:", length(hit), "\n")
}

## ---------- Aggregate ----------
indi_all <- unique(rbindlist(res_indi))       # approximate dedup (same case across quarters)
out <- indi_all[, .(n = .N), by = indi_pt][order(-n)]
fwrite(out, file.path(OUT_DIR, "5_indication_NAION_semaglutide.csv"))
cat("\n===== Indication composition of NAION cases (top 20) =====\n")
print(head(out, 20))

## Diabetes vs weight-management grouping
dm_pt <- c("Type 2 diabetes mellitus", "Diabetes mellitus", "Type 1 diabetes mellitus",
           "Diabetes mellitus inadequate control", "Glucose tolerance impaired",
           "Glycosylated haemoglobin increased", "Blood glucose increased",
           "Insulin resistance")
wt_pt <- c("Weight control", "Obesity", "Weight decreased", "Overweight")
n_dm <- indi_all[indi_pt %chin% dm_pt, uniqueN(primaryid)]
n_wt <- indi_all[indi_pt %chin% wt_pt, uniqueN(primaryid)]
cat("\nAmong NAION cases: diabetes indication", n_dm, "| weight-management indication", n_wt, "\n")

brand_all <- unique(rbindlist(res_brand), by = "primaryid")
cat("\n===== Brand stratification (NAION reports) =====\n")
cat("Ozempic :", sum(brand_all$Ozempic),  "\n")
cat("Wegovy  :", sum(brand_all$Wegovy),   "\n")
cat("Rybelsus:", sum(brand_all$Rybelsus), "\n")
cat("Total (cases):", nrow(brand_all), "\n")