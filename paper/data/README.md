# Laboratory and clinical reference data used for calibration and validation

Every value below was transcribed from the published text or tables of the cited article
(retrieved through PubMed / PubMed Central and ClinicalTrials.gov); the column `source`
gives the location in the article. Units: hours for time, micromolar for concentrations.

| file | what | source |
|---|---|---|
| `iyer2025_u2os_fates.csv` | fates of U2OS cells tracked by time-lapse microscopy for three days after cisplatin at 7, 10 and 13 µM (Low, Medium, High) | Iyer et al. 2025, PLoS Comput Biol 21(9):e1013446, Table 2 and Methods (doi:10.1371/journal.pcbi.1013446) |
| `iyer2025_hct116.csv` | HCT116 lineage experiment at a single cisplatin concentration near the IC50: cells present at drug addition, deaths, and end-fate lineage correlations | same article, Table 1, Fig 1 and Fig 5a |
| `shaffer2017_wm989.csv` | WM989-A6 melanoma: frequency of pre-resistant cells, colony enrichment of sorted EGFR-high cells, drug concentration | Shaffer et al. 2017, Nature 546:431 (doi:10.1038/nature22794) |
| `lasri2020_n15.csv` | N15-0385 patient-derived glioblastoma cells: doubling time and temozolomide concentrations tested | Lasri et al. 2020, R Soc Open Sci 7:191243 (doi:10.1098/rsos.191243) |
| `temozolomide_pk.csv` | temozolomide population pharmacokinetics in glioma patients | Ostermann et al. 2004, Clin Cancer Res 10:3728 (doi:10.1158/1078-0432.CCR-03-0807) |
| `clinical_schedules.csv` | schedules of the RTOG 0525 (NCT00304031) and SWOG S1320 (NCT02196181) trials and their reported outcomes | ClinicalTrials.gov records; Gilbert et al. 2013, J Clin Oncol 31:4085; Algazi et al. 2020, Nat Med 26:1564 |
