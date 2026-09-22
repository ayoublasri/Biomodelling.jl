---
title: "Whether treatment holidays help depends on three measurable properties of resistant cells: a single-cell simulation of two clinical dosing trials"
author:
  - "Ayoub Lasri"
date: "Berache Limited, Dublin, Ireland. Correspondence: alasri@berache.com"
abstract: |
  A melanoma trial and a set of melanoma xenografts gave opposite verdicts on intermittent dosing, and a glioblastoma trial found no benefit from a dose-dense temozolomide regimen designed to deplete MGMT. We ask what a mechanistic single-cell model has to assume to produce each of those outcomes. Resistance here is a heritable expression state with single-cell kinetics rather than a compartment with an assumed rate: a slowly switching promoter whose protein sets a cell's death hazard, inherited through division in cells that grow, replicate their genes and divide. Whether treatment holidays help is decided by three properties of resistant cells — their growth under drug, their fitness without it, and whether the drug itself stabilises the resistant state — and setting those three reproduces both the clinical and the preclinical melanoma outcome, the preclinical one as a lower end-of-simulation burden rather than as a longer time to loss of control. In glioblastoma, modelling MGMT as the suicide enzyme it is creates a dose-dense advantage in unmethylated tumours that the trial did not find, which bounds the per-lesion rate at which temozolomide consumes MGMT; in methylated tumours the model already diverges from the trial with MGMT stable, on cumulative dose alone, so the null there bounds the model's dose response instead. The mechanisms are set in each case rather than inferred from the outcomes, so these are demonstrations that the rankings follow from cell properties, not claims that particular tumours carry them. What the framework adds is that those properties are measurable in lineage-tracing and growth experiments, and can therefore be bounded before a schedule is chosen.
---

## Author summary

Should a cancer drug be given continuously, or with planned breaks? Trials disagree. In melanoma, a randomised trial found continuous dosing better, while mouse experiments had pointed the other way. In glioblastoma, a regimen giving twice as much drug, to exhaust a repair protein that protects cells, did no better than the standard one.

We asked what a model of individual cells has to assume to produce each result. Here a cell is resistant because of which genes it happens to be expressing, a state it passes to its daughters, rather than because it belongs to an abstract "resistant compartment". That makes the assumptions concrete enough to measure.

Three properties decide whether breaks help: how fast resistant cells grow under the drug, how well they do without it, and whether the drug locks the resistant state in place. Set those three, and the model gives either answer, which is why the trials disagreeing is not a paradox. For glioblastoma, the trial's negative result limits how fast the drug uses up the repair protein, but only where that protein is common; where it is rare the model disagrees with the trial for a different reason, and we say so.
