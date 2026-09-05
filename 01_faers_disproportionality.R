## ============================================================================
## Arm 1 — FAERS disproportionality analysis
## GLP-1 receptor agonists x ocular adverse events (main pipeline)
##
## Pipeline: parse quarterly ASCII -> MedDRA standardize -> dedup -> slim cache
##           -> global dedup -> signal metrics (ROR/PRR/IC) -> yearly trends
##           -> time-to-onset -> indication profile
##
## First-time setup:
##   install.packages("BiocManager"); BiocManager::install("faers")
##   install.packages("data.table")
## ============================================================================
library(faers)
library(data.table)

## ---------- 0. Configuration (EDIT THESE PATHS) ----------
FAERS_DIR  <- "/path/to/faers"      # folder containing faers_ascii_YYYYqN subfolders
MEDDRA_DIR <- "/path/to/MedDRA_27_1_English"  # unzipped; point to MedAscii subfolder if an error occurs
CACHE_DIR  <- file.path(FAERS_DIR, "rds_cache")  # per-quarter slim cache (resumable)
OUT_DIR    <- file.path(FAERS_DIR, "output")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

YEARS     <- 2018:2026
QUARTERS  <- c("q1", "q2", "q3", "q4")
ROLE_KEEP <- "PS"   # primary analysis: primary suspect only; sensitivity: c("PS", "SS")

## Drug groups (matched against both drugname and prod_ai, captures fixed combinations)
glp1_members <- c(
  Semaglutide  = "SEMAGLUTIDE|OZEMPIC|WEGOVY|RYBELSUS",
  Tirzepatide  = "TIRZEPATIDE|MOUNJARO|ZEPBOUND",
  Liraglutide  = "LIRAGLUTIDE|VICTOZA|SAXENDA",
  Dulaglutide  = "DULAGLUTIDE|TRULICITY",
  Exenatide    = "EXENATIDE|BYETTA|BYDUREON",
  Lixisenatide = "LIXISENATIDE|ADLYXIN|LYXUMIA"
)
drug_groups <- c(
  list(GLP1RA_all = paste(glp1_members, collapse = "|")),
  as.list(glp1_members),
  SGLT2i       = "EMPAGLIFLOZIN|JARDIANCE|DAPAGLIFLOZIN|FARXIGA|FORXIGA|CANAGLIFLOZIN|INVOKANA|ERTUGLIFLOZIN|STEGLATRO",
  DPP4i        = "SITAGLIPTIN|JANUVIA|SAXAGLIPTIN|ONGLYZA|LINAGLIPTIN|TRADJENTA|ALOGLIPTIN|NESINA",
  Metformin    = "METFORMIN|GLUCOPHAGE",
  Atorvastatin = "ATORVASTATIN|LIPITOR"   # negative-control drug
)

## Ocular outcome PTs (matched against standardized meddra_pt, case-insensitive)
pt_main <- toupper(c(
  "Optic ischaemic neuropathy",                        # primary outcome (NAION PT)
  "Non-arteritic anterior ischaemic optic neuropathy",
  "Anterior ischaemic optic neuropathy",
  "Ischaemic optic neuropathy"
))
pt_secondary <- toupper(c(
  "Blindness", "Visual impairment", "Visual acuity reduced", "Vision blurred",
  "Retinal vein occlusion", "Retinal artery occlusion", "Optic neuropathy",
  "Papilloedema", "Amaurosis", "Amaurosis fugax", "Diabetic retinopathy"
))
pt_all <- unique(c(pt_main, pt_secondary))

## ---------- 1. Per quarter: parse -> standardize -> dedup -> slim cache ----------
load_quarter <- function(y, q) {
  fq <- file.path(FAERS_DIR, sprintf("faers_ascii_%d%s", y, q))
  ## Prefer the unzipped folder; fall back to faers() (recognizes zips in dir,
  ## or downloads from FDA if none exist)
  tryCatch(
    faers_parse(fq, format = "ascii"),
    error = function(e) faers(y, q, dir = FAERS_DIR)
  )
}

for (y in YEARS) for (q in QUARTERS) {
  tag <- sprintf("%d%s", y, q)
  rds <- file.path(CACHE_DIR, paste0(tag, ".rds"))
  if (file.exists(rds)) next                       # resume from checkpoint
  if (!dir.exists(file.path(FAERS_DIR, sprintf("faers_ascii_%d%s", y, q))) &&
      !file.exists(file.path(FAERS_DIR, sprintf("faers_ascii_%d%s.zip", y, q)))) next
  message(">>> Processing ", tag)
  ok <- tryCatch({
    x <- load_quarter(y, q)
    x <- faers_standardize(x, MEDDRA_DIR)          # MedDRA PT standardization
    x <- faers_dedup(x)                            # within-package dedup
    slim <- list(
      demo = faers_get(x, "demo")[, .(primaryid, caseid, caseversion, fda_dt, event_dt)],
      drug = faers_get(x, "drug")[, .(primaryid, drug_seq, role_cod, drugname, prod_ai)],
      reac = faers_get(x, "reac")[, .(primaryid, pt = toupper(meddra_pt))],
      indi = faers_get(x, "indi")[, .(primaryid, indi_pt)],
      ther = tryCatch(faers_get(x, "ther")[, .(primaryid, dsg_drug_seq, start_dt)],
                      error = function(e) NULL)
    )
    saveRDS(slim, rds)
    rm(x, slim); gc()
    TRUE
  }, error = function(e) { message("!! ", tag, " failed: ", conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) next
}

## ---------- 2. Merge + second-pass global dedup across quarters ----------
files <- list.files(CACHE_DIR, "\\.rds$", full.names = TRUE)
slim <- lapply(files, readRDS)
demo <- rbindlist(lapply(slim, `[[`, "demo"))
drug <- rbindlist(lapply(slim, `[[`, "drug"))
reac <- rbindlist(lapply(slim, `[[`, "reac"))
indi <- rbindlist(lapply(slim, `[[`, "indi"))
ther <- rbindlist(lapply(slim, function(x) x$ther), fill = TRUE)
rm(slim); gc()

demo[, caseversion := as.integer(caseversion)]
setorder(demo, caseid, -caseversion, -fda_dt)
keep_id <- demo[!duplicated(caseid), primaryid]    # latest version of each case
demo <- demo[!duplicated(caseid)]
drug <- drug[primaryid %chin% keep_id]
reac <- reac[primaryid %chin% keep_id]
indi <- indi[primaryid %chin% keep_id]
ther <- ther[primaryid %chin% keep_id]
demo[, year := as.integer(substr(fda_dt, 1, 4))]
N_total <- nrow(demo)

## ---------- 3. Drug-case table ----------
drug[, namekey := toupper(paste(drugname, prod_ai))]
drug <- drug[role_cod %chin% ROLE_KEEP]
dl <- rbindlist(lapply(names(drug_groups), function(g) {
  sub <- drug[grepl(drug_groups[[g]], namekey)]
  if (!nrow(sub)) return(NULL)
  unique(sub[, .(primaryid, group = g)])
}))

## ---------- 4. Signal metrics (ROR / PRR / IC) ----------
N_total <- uniqueN(demo$primaryid)   
reac_u <- unique(reac[pt %chin% pt_all, .(primaryid, pt)])
n_g  <- dl[, .(n_group = uniqueN(primaryid)), by = group]
n_pt <- reac_u[, .(n_pt = uniqueN(primaryid)), by = pt]
a_tab <- merge(reac_u, dl, by = "primaryid", allow.cartesian = TRUE)[
  , .(a = uniqueN(primaryid)), by = .(group, pt)]
sig <- merge(CJ(group = names(drug_groups), pt = pt_all), a_tab,
             by = c("group", "pt"), all.x = TRUE)[is.na(a), a := 0]
sig <- merge(sig, n_g,  by = "group", all.x = TRUE)
sig <- merge(sig, n_pt, by = "pt",    all.x = TRUE)
sig[is.na(n_pt), n_pt := 0]
sig[, `:=`(b = n_group - a, c = n_pt - a)]
sig[, d := N_total - a - b - c]
sig[, ROR    := ((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5))]
sig[, ROR_lo := exp(log(ROR) - 1.96 * sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5)))]
sig[, ROR_hi := exp(log(ROR) + 1.96 * sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5)))]
sig[, PRR    := ((a + 0.5) / (a + b + 1)) / ((c + 0.5) / (c + d + 1))]
sig[, E      := (n_group * n_pt) / N_total]
sig[, IC     := log2((a + 0.5) / (E + 0.5))]
sig[, IC025  := IC - 1.96 * sqrt((1/log(2))^2 * (1/(a+0.5) + 1/(E+0.5)))]
setorder(sig, group, -a)
fwrite(sig, file.path(OUT_DIR, "1_signal_table.csv"))

## ---------- 5. Yearly trends ----------
dl_y <- merge(dl, demo[, .(primaryid, year)], by = "primaryid")
rx_y <- merge(reac_u, demo[, .(primaryid, year)], by = "primaryid")
yr <- merge(rx_y, dl_y[, .(primaryid, group)], by = "primaryid", allow.cartesian = TRUE)
fwrite(yr[, .(a = uniqueN(primaryid)), by = .(year, group, pt)],
       file.path(OUT_DIR, "2_yearly_group_pt.csv"))       # year x drug group x PT
fwrite(dl_y[, .(n_group = uniqueN(primaryid)), by = .(year, group)],
       file.path(OUT_DIR, "2_yearly_group_total.csv"))    # year x drug group totals
fwrite(rx_y[, .(n_pt = uniqueN(primaryid)), by = .(year, pt)],
       file.path(OUT_DIR, "2_yearly_pt_total.csv"))       # year x PT totals
fwrite(demo[, .(N = .N), by = year],
       file.path(OUT_DIR, "2_yearly_total.csv"))          # whole-database yearly totals

## ---------- 6. Time-to-onset (GLP-1 RA x primary ocular outcome) ----------
eye_main  <- unique(reac[pt %chin% pt_main, .(primaryid, pt)])
glp1_drug <- drug[grepl(paste(glp1_members, collapse = "|"), namekey),
                  .(primaryid, drug_seq, namekey)]
tto <- merge(glp1_drug, ther, by.x = c("primaryid", "drug_seq"),
             by.y = c("primaryid", "dsg_drug_seq"), all.x = TRUE)
tto <- merge(tto, eye_main, by = "primaryid")
tto <- merge(tto, demo[, .(primaryid, event_dt)], by = "primaryid")
tto <- tto[nchar(start_dt) == 8 & nchar(event_dt) == 8]
if (nrow(tto)) {
  tto[, tto_days := as.integer(as.Date(event_dt, "%Y%m%d") - as.Date(start_dt, "%Y%m%d"))]
  for (g in names(glp1_members)) tto[, (g) := grepl(glp1_members[[g]], namekey)]
  fwrite(tto[, .(primaryid, pt, tto_days, Semaglutide, Tirzepatide, Liraglutide,
                 Dulaglutide, Exenatide, Lixisenatide)],
         file.path(OUT_DIR, "3_tto_glp1_eye.csv"))
}

## ---------- 7. Indication profile of semaglutide reports ----------
sema_ids <- dl[group == "Semaglutide", primaryid]
fwrite(indi[primaryid %chin% sema_ids, .(n = .N), by = indi_pt][order(-n)][1:50],
       file.path(OUT_DIR, "4_indication_semaglutide.csv"))

message("Done! Output directory: ", OUT_DIR)

