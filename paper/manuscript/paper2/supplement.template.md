---
title: "Supplementary Information"
author: "Ayoub Lasri"
date: "Berache Limited, Dublin, Ireland"
---

**Whether treatment holidays help depends on three measurable properties of resistant cells: a single-cell simulation of two clinical dosing trials.** Supplementary Notes 1–4, Supplementary Figure 1 and Supplementary Tables 1–2.

# Supplementary Note 1: Cell-cycle-dependent killing in the melanoma schedules

The death hazard used here has no cell-cycle dependence. The companion Article shows that the cisplatin fate data cannot decide whether cells are spared because they inherited a protective state or because they were outside the replication window, so the same question applies to the schedule comparisons: do their rankings rest on that choice? We repeated the melanoma study with three quarters of the death hazard confined to a window around mid-cycle ($\beta = 0.25$, $\varphi_0 = 0.5$, $w = 0.15$), leaving everything else unchanged, and again with the maximal hazard raised so that the cycle-averaged hazard matches the cycle-blind one.

The window position is a modelling choice rather than a claim about the agent. The melanoma study is built on a BRAF/MEK inhibitor, which does not create replication-coupled lesions as a platinum drug does, although it does act on cells traversing the cycle. What these runs test is robustness to a strong cell-cycle gate of any kind, not the phase specificity of a particular drug.

![](paper2_figS1_cycle.png)

**Supplementary Fig. 1 | The melanoma schedules under cell-cycle-gated killing.** **a**, Weeks after randomisation to loss of control for the three mechanisms and the three clinical schedules (mean ± s.d. of four seeds, censored at the end of follow-up), cycle-blind (solid), cycle-gated (mid) and gated at matched hazard (pale).

{{cycle_mel_note}} {{cycle_matched_mel_note}}

# Supplementary Note 2: Sensitivity of the melanoma schedules to the memory of the resistant state

{{memory_scan_note}} {{memory_scan_reconcile}}

# Supplementary Note 3: Parameters of the simulations

Supplementary Table 1. Parameters of every simulation in this Article. The shared model layers are described in the Methods of the companion Article.

| Figure | Model | Parameters |
|---|---|---|
| 1 | melanoma-like | $k_{\mathrm{tx}} = 3$, $k_{\mathrm{dm}} = 0.1$, $k_{\mathrm{tl}} = 0.4$, $k_{\mathrm{dp}} = 0.02$ per h; net doubling 4 weeks, memory 5 net doublings, pre-resistant fraction 0.005 (1:200, above the traced 1:1,000 to 1:10,000 of @emert2021; Methods), $h_{\max} = 0.0023$ per h, $\mathrm{EC}_{50} = 0.3$, $m = 2$, protection $K = 150$, $q = 4$, growth arrest $\mathrm{IC}_{50} = 0.3$, $m_g = 2$, with the same protection, $k_{\mathrm{off}} \to k_{\mathrm{off}}/(1 + 9d)$, fitness cost 0.5, partial protection: additional unprotected growth inhibition with $\mathrm{IC}_{50} = 1$; 2,000 founders with promoter states drawn from the stationary distribution, 60 weeks, dt 4 h; progression from the nadir at 1.73 × the running minimum after the 8-week lead-in, loss of control at 1.2 × the pre-treatment size; cycle-gated variant (Supplementary Fig. 1) adds $\beta = 0.25$, $\varphi_0 = 0.5$, $w = 0.15$, and is repeated with $h_{\max} = 0.0045$ per h so that the cycle-averaged hazard matches the cycle-blind one; the memory of the resistant state is scanned over 2, 3.5, 5, 8 and 12 net doublings with three seeds per point |
| 2 | MGMT model | same expression kinetics; net doubling 40 days, memory 4 net doublings, MGMT-expressing fraction 0.01 or 0.30, $h_{\max} = 0.03$ per h at peak, $\mathrm{EC}_{50} = 0.4$ of the standard bolus peak, $m = 2$, $K = 150$, $q = 4$; elimination half-life 2.1 h; stoichiometric MGMT consumption at rate $k\,d\,V\,c/(c + K_m)$ with $k = 300$, $K_m = 150$, $k$ chosen so that the consumption rate at the standard bolus peak matches the first-order parameterisation it replaces rather than fitted to data, and with no lesion pool as a state variable, so lesions are assumed to form in proportion to dose and to be repaired at once; 1,500 founders at the stationary promoter distribution, six 28-day cycles, dt 1 h |

# Supplementary Note 4: Clinical and laboratory reference data

Supplementary Table 2. Reference values used by the case studies (transcribed from the cited articles; files in `paper/data/`).

| Quantity | Value | Source |
|---|---|---|
| Pre-resistant melanoma cells traced back from resistant fates | initial frequency ~1:1,000 to 1:10,000 (context, not a value used: the melanoma study here uses 1:200; Methods) | @emert2021 |
| Pre-resistant melanoma cells | 1:50 to 1:500 per marker; EGFR-high cells give 7.9 ± 0.9 fold more resistant colonies | @shaffer2017 |
| N15-0385 glioblastoma doubling time | 50 h | @lasri2020 |
| Temozolomide elimination half-life in plasma | 1.8 h (single dose); 2.1 h (population model of plasma and cerebrospinal fluid; the value used here) | @rudek2004; @ostermann2004 |
| Temozolomide penetration of the cerebrospinal fluid | exposure 20% of plasma exposure | @ostermann2004 |
| MGMT / alkyltransferase depletion in peripheral blood mononuclear cells | $-63\%$ at 14 days, $-73\%$ at 21 days on protracted schedules; nadir 18.0 ± 2.26% of initial on a compressed 1,000 mg/m² schedule | @tolcher2003; @middleton2000 |
| Tumour MGMT activity in orthotopic GBM43 xenografts | depleted by day 6 on both schedules; still suppressed at day 22 only on the 21-day schedule; back to baseline in both by day 29 | @robinson2010 |
| RTOG 0525 regimens | 150-200 mg/m² days 1-5 vs 75-100 mg/m² days 1-21 of 28-day cycles; median OS 16.6 vs 14.9 months | @gilbert2013 |
| SWOG S1320 regimens | continuous vs 3 weeks off / 5 weeks on after an 8-week lead-in; median PFS 9.0 vs 5.5 months (HR 1.36 intermittent:continuous, $P = 0.063$; the trial pre-specified two-sided $\alpha = 0.2$ and 80% confidence intervals); median OS 29.2 months in both arms, a secondary end point the trial was not powered for | @algazi2020 |

