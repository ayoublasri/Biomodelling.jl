---
title: "Mechanistic simulation of heritable expression states, cell division and drug response in single-cell populations"
author:
  - "Ayoub Lasri"
date: "Independent researcher, Dublin, Ireland. Correspondence: alasri@berache.com"
abstract: |
  Heritable differences in gene expression decide which cancer cells survive a drug, yet the simulators used to benchmark single-cell methods model neither growth nor division, and stochastic simulators of dividing cells carry no drug layer. Here we present Biomodelling.jl 2.0, in which stochastic reaction kinetics run inside growing, dividing, drug-treated cells. Exact, hybrid and adaptive kernels couple to exponential growth, volume-scaled transcription, gene replication, division with partitioning, dose schedules with state-dependent killing, a lineage record, observation models and likelihood-free inference. Expression memory, cell-size scaling and cell-cycle-dependent bursting emerge from this physiology rather than imposed. One resistance gene with slow promoter switching reproduces the lineage and memory signatures of drug tolerance, and supplies ground truth that benchmarks of memory genes, network inference and perturbation prediction lack. Calibrated to time-lapse measurements of cisplatin-treated cells, it reproduces the dose-invariance of single-cell timing and the kin fate correlations, though six parameters fitted to six fate fractions leave that memory undetermined, and a cycle-gated hazard fits as well. Its schedule optimiser reproduces the opposite rankings of intermittent dosing reported in a melanoma trial and in xenografts, and shows the ranking to turn on two measurable cell properties: growth under drug and fitness without it.
---

## Author summary

Two cancer cells with the same genome, given the same drug, can meet opposite fates: one dies and the other survives. Much of that difference comes from how strongly each cell happens to be expressing particular genes when treatment begins, and cells pass those expression states to their daughters when they divide. The simulators used to test single-cell analysis methods do not model growth or division, so they cannot produce such inherited states; the simulators that do model dividing cells carry no drug.

We built Biomodelling.jl 2.0, in which chemical reactions play out inside cells that grow, copy their genes, divide, and meet drugs on realistic schedules. Expression memory, cell-size effects and cell-cycle dependence emerge from that physiology rather than being assumed. We checked the simulator against cases whose answers are known exactly in mathematics, fitted it to published time-lapse recordings of cells under chemotherapy, and asked it to predict a concentration it had not been fitted to.

We then put a clinical question to it: should a drug be given continuously, or with breaks? The answer turns on two properties of resistant cells that experiments can measure.
