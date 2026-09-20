# LinkedIn post — video launch

**Attach:** `biomodelling-showcase.mp4` (1080×1080, 60 s, silent — designed to read without sound)

---

Most single-cell simulators can't let a cell grow or divide.

That sounds like an implementation detail. It isn't. Heritable, non-genetic differences in
gene expression are a large part of why one cancer cell survives a drug and the cell beside it
doesn't — and those states come from growth, gene replication, division, and what daughters
inherit. A simulator without that physiology can't generate the phenomenon the methods it is
used to benchmark were built to detect.

Biomodelling.jl 2.0 runs stochastic reaction kinetics inside growing, dividing, drug-treated cells.

Things that fall out of the physiology instead of being assumed:
— expression memory across generations
— cell-size scaling of transcript counts
— cell-cycle-dependent bursting
— persister-like fractional killing, with sisters sharing their lineage's fate 4.8× more often
  than chance

It is calibrated, not just illustrative: fitted to published time-lapse measurements of
cisplatin-treated U2OS cells, with one concentration held out.

And it reports what the data cannot determine. The fate counts pin the death parameters to
about twofold but leave the promoter switching rates unidentified — memories from one to nine
generations, paired with resistant fractions from 1% to 25%, all fit the same measurements
within noise. Refitting with cell-cycle-dependent killing fits just as well, so those counts
cannot separate inherited protection from where a cell sat in its cycle when the drug arrived.

As a schedule testbed it reproduces the opposite intermittent-dosing outcomes of a melanoma
trial and of xenografts, once you set how resistant cells grow with and without drug. And
modelling MGMT as the suicide enzyme it actually is predicts a dose-dense temozolomide
advantage that RTOG 0525 did not find — which bounds how fast the drug can really be consuming
MGMT in tumours.

Open source, MIT licensed, white paper included.

github.com/ayoublasri/Biomodelling.jl

#computationalbiology #singlecell #drugdiscovery #julialang #cancerresearch #systemsbiology
