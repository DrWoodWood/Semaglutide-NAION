# Semaglutide and NAION: Pharmacovigilance Signal or Notoriety Bias?

Code for the two-arm study:

- **Arm 1 — Pharmacovigilance:** Disproportionality analysis of NAION reports for
  semaglutide in FDA FAERS (2018–2026), including yearly reporting trends,
  notoriety-bias assessment, time-to-onset analysis, and indication/brand
  composition of NAION cases.
- **Arm 2 — Genetics:** Drug-target Mendelian randomization of *GLP1R*
  expression (eQTLGen whole-blood cis-eQTLs) against glaucoma-spectrum ocular
  traits used as NAION proxies, with HbA1c as positive control and skin tanning
  as negative control.

## Repository structure

```
R/
  01_faers_disproportionality.R   Arm 1 main pipeline: parse FAERS quarterly ASCII,
                                  dedup, ROR/PRR/IC025, yearly trends, TTO
  02_naion_indication_brand.R     Arm 1 supplement: NAION indication composition
                                  (diabetes vs weight loss) and brand split
                                  (Ozempic/Wegovy/Rybelsus)
  03_glp1r_instruments.R          Arm 2 instrument selection: GLP1R cis-eQTLs +
                                  LD clumping -> glp1r_lead_snps.csv
  04_drug_target_mr.R             Arm 2 TwoSampleMR pipeline: harmonisation,
                                  IVW/Egger/weighted-median, heterogeneity,
                                  pleiotropy, leave-one-out, sensitivity
  05_figures.R                    Publication figures: Figure 1 (yearly NAION
                                  reporting trend) and Figure 2 (two-arm forest
                                  plot); reads results CSVs, falls back to the
                                  manuscript's key values for standalone use
results/                          Output CSVs (created by the scripts)
figures/                          Output figures (created by 05_figures.R)
```

Run the scripts in order `01 → 04`. Scripts 01–02 are independent of 03–04.

## Requirements

R (>= 4.2) with packages:

- `data.table`
- [`faers`](https://cran.r-project.org/package=faers) (Bioconductor-style FAERS parsing)
- `ieugwasr` (LD clumping and OpenGWAS access)
- `TwoSampleMR`
- `ggplot2` (+ optional `patchwork` to combine Figure 2 panels)

## Data sources

| Data | Source | Notes |
|---|---|---|
| FAERS quarterly ASCII files | [FDA FAERS Quarterly Data](https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files) | Download and unzip all quarters 2018Q1–2026Qx into one directory; set `FAERS_DIR` in scripts 01/02 |
| MedDRA dictionary | MSSO (licensed) | Requires a MedDRA license; set `MEDDRA_DIR` |
| eQTLGen significant cis-eQTLs | [eqtlgen.org/cis-eqtls.html](https://eqtlgen.org/cis-eqtls.html) | Download the *Significant cis-eQTLs* file (not the full summary statistics); set `EQTL_FILE` in script 03 |
| Outcome GWAS | OpenGWAS / IEU | Accessed by dataset ID; requires an OpenGWAS JWT token: `Sys.setenv(OPENGWAS_JWT = "...")` |

## Outputs

Arm 1 (scripts 01–02):

- `1_signal_metrics.csv` — ROR/PRR/IC for NAION across drug groups
- `2_yearly_cases.csv`, `3_yearly_reports.csv`, `4_yearly_rate.csv` — yearly trends
- `5_indication_NAION_semaglutide.csv` — indication and brand composition
- TTO summary (Weibull shape parameter)

Arm 2 (scripts 03–04):

- `glp1r_lead_snps.csv` — clumped GLP1R instruments
- `B1_mr_results.csv` — IVW / Egger / weighted median / Wald ratio estimates
- `B1_heterogeneity.csv`, `B1_pleiotropy.csv`, `B1_leaveoneout.csv`
- `B1_mr_results_noLowN.csv` — sensitivity excluding rs114977861 (low-N instrument)

## Notes on reproducibility

- FAERS data are periodically updated by the FDA; re-running the pipeline on a
  later data release may yield slightly different counts.
- OpenGWAS datasets occasionally move or are re-harmonised; if an outcome ID
  returns no SNPs, check `gwasinfo()` for the current status of the dataset.
- All analyses use global deduplication by `caseid` + max(`caseversion`) and a
  0.5 continuity correction for disproportionality metrics.

## License

This project is released under the [MIT License](LICENSE).

## Citation

If you use this code, please cite the accompanying manuscript
(*Semaglutide and NAION: Pharmacovigilance Signal or Notoriety Bias?*, under review).
