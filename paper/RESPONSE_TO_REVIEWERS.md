# Response to the reviewer

We thank the reviewer for a detailed and unusually specific report. Every enumerated
point is addressed below, with the change made and where to find it. Three of the
points changed conclusions rather than wording, and we say so plainly where that
happened.

The revision also adds two capabilities to the software that the review made
necessary, a mother-machine population control and general interdivision-time
distributions, and a validation of the population layer against exact solutions
that now runs as part of the test suite on every commit.

---

## Required (blocking)

### 1. Reconcile all cisplatin error metrics and every text-versus-figure mismatch

**Issue 1, the three error metrics.** The reviewer is right that a penalised
objective cannot be smaller than the error it contains. The three numbers were
three different quantities, computed on different seed sets, none of which said so.
They are now defined once and labelled everywhere:

| Quantity | Definition | Value |
|---|---|---|
| Training error | root-mean-square residual of the seed-averaged fate fractions, pooled over all six fractions of the two training concentrations, four seeds | 0.077 |
| Held-out error | the same over the three fractions at 10 µM | 0.093 |
| Calibration objective (colour scale of Fig. 7f) | what the optimiser minimised: single-seed training error plus the memory prior | 0.067 = 0.061 + 0.006 |
| Model comparison (Supplementary Note 12) | the single-seed training error, used identically for both killing mechanisms | 0.061 against 0.064 |

The Results now state all four and say why they differ. The pooled definition
replaces the previous mean-of-two-per-concentration errors, which is why the
training figure moves from 0.071 to 0.077.

**Issue 2, Fig. 3f.** The claim of "the expected factor of two in burst frequency"
is withdrawn. What the panel shows is now reported from the panel's own numbers:
the fitted transcription rate rises 2.0-fold between the one- and two-copy subsets,
of which 1.4-fold is the larger volume of cells that have replicated, leaving
1.4-fold per unit volume rather than twofold; the fitted activation rate also rises
1.4-fold. We say that part of the shortfall is misspecification, since two
independently switching alleles are not one telegraph gene, and give the fitted
values against the truth (0.5, 1.5 and 2 per unit volume per allele against 0.59,
1.89 and 2.41 for the one-copy subset). The panel is now presented as a
demonstration that a cycle-blind fit is biased, not as a measurement of the
copy-number effect. The legend also states the normalisation that made the
transcription axis read 2.4 and 3.4.

**Issue 3, Fig. 5a.** The text reported the mean of the slowest group, 29.6, next to
a plotted maximum of 44. Both are now given: a mean of 29.6 for genes with memory
beyond a hundred cell-cycle times, with the highest single gene at 44.5 (38.4 on
sequenced counts).

**Issue 7, Fig. 8a.** Corrected, and the reviewer's reading was the right one.
Without a fitness cost both interrupted schedules do keep control longer than
continuous dosing (23.4, 26.1 and 26.7 weeks). The sentence claiming that holidays
"only let sensitive cells regrow and be killed again" is replaced by the actual
mechanism, which was previously only in the Methods: the drug stabilises the
resistant state tenfold by suppressing the promoter inactivation rate, so continuous
exposure keeps cells resistant that a holiday lets revert. The 65% cumulative dose
is now stated as being over the full 60-week simulation.

**Issue 9, Fig. 7c.** The third-cousin correlation is now reported with a confidence
interval over the three seeds, and the text states that the source reports these
correlations graphically and gives no numeric coefficients, so the comparison is of
pattern rather than of value.

### 2. The Harmange 2023 contradiction and the PI3K-versus-DOT1L attribution

Both corrected, and the Fig. 4 narrative now treats the discrepancy directly rather
than implying agreement.

The mechanism attribution was wrong: Harmange et al. identify PI3K and TGF-β
signalling as switching modulators, and DOT1L comes from the Torre et al. CRISPR
screen (Nat. Genet. 2021, 53:76-85). That screen is also the experiment that tests
what is simulated here, and it agrees with the model: inhibiting DOT1L *before*
adding a BRAF inhibitor produced **more** therapy resistance than giving the two
together.

The two experiments differ in what the pretreatment does, and the text now says so.
The perturbation in the model is symmetric: it multiplies the activation and the
inactivation rate by the same factor, so it accelerates relaxation towards the
stationary resistant fraction without moving that fraction, and before the drug its
only effect is to spread the resistant state over more lineages. A modulator that
shifts the balance between the states is a different perturbation: a pretreatment
that drives cells out of the primed state lowers the number of protected cells when
the drug arrives and reduces resistance whether or not it is maintained. The
prediction is therefore stated as specific to perturbations that speed switching
without biasing it, and the text names the measurement that distinguishes the two
cases.

### 3. The temozolomide/MGMT "not measured" claim

Softened as requested and made more specific. The Results now say that the quantity
the argument rests on is the per-lesion molecular consumption rate constant, which
has not been measured directly, and cite the bulk depletion time-courses that have:
−63% at 14 days and −73% at 21 days in peripheral blood mononuclear cells on
protracted schedules (Tolcher et al. 2003), a nadir of 18.0 ± 2.26% of initial levels
on a compressed schedule with O⁶-methylguanine measured in tumour biopsies as well as
blood (Middleton et al. 2000), and tumour MGMT activity in orthotopic GBM43
xenografts depleted by day 6 on both schedules and still suppressed at day 22 only on
the 21-day schedule (Robinson et al. 2010, Br. J. Cancer 103:498-504).

The xenograft result also sharpens the disagreement rather than merely restating the
trial, and the Results now make that point: the deeper and longer tumour depletion
that the dose-dense schedule is supposed to deliver was confirmed in that model, and
survival still did not improve. Either bulk depletion is not the per-cell depletion
that matters for killing, or the benefit of depleting MGMT is offset by something the
model leaves out. All three references are added to the reference data table with
their numbers.

### 4. Table 1, Supplementary Note 10 and the novelty claim

**Table 1 corrections.** SERGIO is Python and returns expression matrices as arrays,
and scMultiSim is R and returns SingleCellExperiment-style objects, so the AnnData
entry for both is changed from ✓ to –. scMultiSim carries a differentiation tree, so
its lineage entry is changed from – to ○, and the caption says cells are related by
a trajectory even though per-cell division and fate are not simulated.

**AgentBasedModeling.jl.** Four entries were understated and are corrected: gene
regulatory networks with known ground truth, heritable promoter states and parameter
inference go from ○ or – to ✓, and gene replication goes from – to ○. The inference
entry matters most: the companion work (Piho & Thomas, Sci. Adv. 2024, 10:eadl4895)
infers kinetic and selection parameters from mother-machine and lineage-tree data by
finite state projection, which the previous "–" denied. That reference is added. The
caption now justifies the remaining ○ entries: gene replication is expressible in the
agent formalism but not provided, and the lineage output traces a cell back to its
ancestor rather than recording a tree with per-cell fates.

**Version 1 overlap.** The Discussion now states precisely what v2 adds over v1:
gene replication and cell-cycle copy number, sizer and adder division rules in place
of a timer, inherited promoter states, the lineage table, the drug and
pharmacokinetic layer, likelihood-free inference and the schedule optimiser.
Supplementary Note 10 gains a row saying that the synthetic-scRNA-seq benchmarking of
Fig. 5 is something v1 already did, and that nothing in that figure which v1 could
also produce is presented as new.

**Novelty reframed.** The Discussion's "three things are new" is replaced by the
reviewer's formulation: no individual mechanism is new, and this is the first package
to combine these mechanisms with joint inference and schedule optimisation in a single
forward model. The Introduction's "none of them can" statements were audited against
the wider list and softened accordingly, and Splatter, SymSim, scDesign3, PhysiCell
and PhysiBoSS are now cited and positioned.

### 5. Validation against the exact solutions, in both lineage and snapshot modes

Done, and it required two additions to the software.

The package gains a `MotherMachine` population control, which keeps one daughter at
random at every division so that N founders give N statistically independent
single-lineage traces, and `AgeTimer` gains gamma-distributed interdivision times,
which covers the Erlang family of the exactly solvable models and, at a coefficient
of variation of one, the memoryless timer. `simulate_population` also accepts the
founders' ages.

Four analytically solvable members of the class were simulated in both modes:
constitutive production, production in bursts of fixed size and a two-state promoter,
each with a memoryless interdivision time; and volume-scaled synthesis from a gene
that replicates at mid-cycle, with a deterministic cycle and exponential volume
growth. For the first three the stationary law is available in closed form in both
settings, with means k/(γ + λ/2) along a lineage and k/(γ + λ) in a population, and
the full distributions follow from the stationary state of one update step. For the
fourth the count distribution at every cycle phase is exactly Poisson and the two
modes differ only in how phases are weighted.

The results are in Supplementary Note 14, Supplementary Fig. 6 and Supplementary
Table 5. Every simulated sample matches the exact law of its own mode below the 99%
critical value for its sample size, and scored against the law of the *other* mode the
same samples give distances an order of magnitude larger, so the comparison has ample
power to tell the two settings apart. The Methods describe the protocol, and the
closed-form lineage and population means and variances are now checked on every run
of the test suite (`test/test_exact_population.jl`).

The Jia & Grima venue is corrected to iScience (2023, 26:105746), and Beentjes et al.
(Phys. Rev. E 2020, 101:032403), Thomas (J. R. Soc. Interface 2017) and Thomas &
Shahrezaei (J. R. Soc. Interface 2021) are cited for the snapshot-versus-lineage
framing.

One incidental finding is worth recording. With an update step that does not divide
the interdivision time exactly in binary, the test `age >= T` falls on a
floating-point boundary and about half the cells divide one step late. It does not
affect the count distributions, but it smears the cycle-phase lattice, and the
validation now uses a step of 1/512 of the interdivision time.

### 6. The cycle-gate comparison at matched cycle-averaged mean hazard

Done, and the reviewer's arithmetic is confirmed. Averaged over a uniform cycle phase
the multiplier is 0.532. Under the phase density these simulations actually realise,
f(φ) = 2/(1+φ)² for a sizer with exponential growth in a growing population, it is
0.508, and we use the latter. A second arm of every cycle-gated comparison now raises
the maximal hazard by 1/0.508 so that the cycle-averaged hazard matches the
cycle-blind one, and both arms are reported. Supplementary Note 13 states the two
averages explicitly and says that half of the shallower dose response under the gate
is simply less killing.

The conclusion the comparison was drawn for does not depend on the matching: at
matched mean hazard the melanoma schedule ranking is essentially the cycle-blind one
again, and under partial protection continuous dosing still outlasts the intermittent
schedule, which is the trial-consistent ordering.

### 7. Seed replication on the schedule scans, and Supplementary Table 1

Both done. Every point of the Fig. 4h schedule scan and of its cycle-gated
counterpart is now the mean of five independent seeds, and the text states the gap
between the best schedule and the next best against the standard error on that
difference, so that where the scan identifies a region of good schedules rather than
a single best one it says so.

The melanoma memory of five net doublings (issue 13) now has the sensitivity scan the
reviewer asked for, over two to twelve net doublings with three seeds at each point
(Supplementary Note 15). The ranking is unchanged throughout; what changes is how
much the holidays are worth, and the advantage of the intermittent schedule over
continuous dosing falls from about nine weeks at two doublings to under two weeks at
twelve. The assumed value therefore sits at the conservative end.

Supplementary Table 1 is restored, as a table of the correctness tests of the kernels
and of the population layer with their enforced tolerances, including the new
exact-law tests.

### 8. The limitation about second-cousin correlations

Added, prominently, in both the Fig. 7 narrative and the Discussion. Cells here draw
their interdivision times independently, so no cycle timing is inherited and every
fate correlation between related cells is carried by the inherited promoter state.
Cell-cycle duration is in fact heritable over several generations and its inheritance
alone produces strong cousin and second-cousin correlations even where mother-daughter
correlations are weak (Sandler et al. 2015; Kuchen et al. 2020; Chakrabarti et al.
2018, all now cited). Combined with a cycle-gated hazard that is a complete competing
explanation for exactly the correlations the calibration reproduces, and the model as
built cannot generate it, so the fit has no choice but to load all of that correlation
onto expression memory. The Discussion names heritable cycle timing as the extension
this points to.

---

## Recommended (non-blocking)

**9. Corigliano et al. 2025 compared directly against Fig. 8.** Done, in the Results
where the optimiser output is reported. The agreement is partial and the difference is
informative: where tolerance is drug-induced the optimiser here also selects an
intermediate dose over the highest one, for the same reason, but release periods help
only when the drug changes the kinetics of the resistant state, either by stabilising
it or by inducing it at a saturating rate. With a state whose decay rate is the same
with and without drug, a holiday only postpones the death of reverting cells.

**10. SWOG S1320 post-progression survival and the drug-addiction mechanism.** Both
added to the Discussion: the progression-free survival advantage of continuous dosing
did not carry through to survival after progression, which favoured the intermittent
arm; and the resistant xenograft tumours were drug-addicted, regressing on withdrawal
rather than merely growing more slowly, which is the extreme of the fitness cost the
second mechanism represents.

**11. The Fig. 5b AUPR inversion.** Explained where it occurs. Dividing counts by the
library size removes cell-to-cell scale, and in a dividing population that scale is
dominated by volume and gene copy number, so the normalisation acts as a coarse
correction for the growth confound that dividing by the true volume does not make. The
text notes that both the fixed-volume ceiling and the cell-cycle regression are above
it, so this is not evidence that sequencing helps.

**12. Runtime scaling qualified per kernel.** Done. The claim of linear scaling in the
number of reactions is now made only for the hybrid kernel, with both series quoted,
and the reason given: the direct method draws one event at a time while the hybrid
kernel leaps over the non-critical channels in a fixed number of operations.

**13. The abstract aligned to the identifiability analysis.** Rewritten. It now says
the model reproduces the dose-invariance of single-cell timing and the kin fate
correlations, that six parameters fitted to six fate fractions leave the memory behind
them undetermined, and that a cycle-gated hazard fits as well. It is still 200 words.

---

## Other enumerated points

**Issue 11, the lesion pool and the rate constant.** Both stated, in the Results where
the mechanism is introduced and again in the Discussion. Lesions are not a state
variable: they are taken to form in proportion to dose and to be repaired at once, so
"one molecule per lesion" enters only through the rate of consumption and no
unrepaired damage accumulates. The rate constant is not fitted to data; it is set so
that consumption at the standard bolus peak matches the first-order parameterisation
it replaces.

**Issue 12, the 30% MGMT-expressing fraction.** Flagged where the prediction is made,
not only in the limitations: the two fractions are illustrative values that bracket
low and high expression, promoter methylation shifts expression as a graded quantity
rather than switching a discrete subpopulation on, and the whole dose-dense prediction
rests on treating it as two states.

**Issue 14, the small-network caveat.** Now adjacent to each result rather than only in
the Discussion, in both the network-inference and the perturbation-prediction
sentences.

**Iyer et al. wording.** The correlations are described as first- and second-cousin
throughout, and the text and the reference table both state that the source gives them
graphically with no numeric coefficients.

**Temozolomide half-life.** Both values are given: 1.8 h in plasma (Rudek et al. 2004)
and 2.1 h in the cerebrospinal-fluid population model used here (Ostermann et al.
2004), in the Results, the Fig. 8 legend and the reference table.

**Shaffer Rewind frequencies.** Added where the barcode signature is discussed: the
follow-up method traces resistant fates back to drug-naive precursors at an initial
frequency of about 1 in 1,000 to 1 in 10,000 cells (Emert et al. 2021), the regime the
melanoma case study is set to.

**AI-authorship disclosure.** Moved to a dedicated Methods statement, "Use of
artificial intelligence", and reflected in the author contributions. It specifies what
the model did (the Julia source, the simulation and plotting scripts, the test suite
and the first draft of the text, from the author's specification) and what the author
verified independently (every scientific question, the choice of models, parameters
and reference data, all numerical results against the generated output tables, and
every claim against its cited source). A second paragraph states that no display item
is generated by a generative model: every figure and table is produced
deterministically from simulation output by the scripts in `paper/`, and every number
in the Results and Supplementary Information is substituted from the output tables by
`paper/manuscript/fill_numbers.py`.

---

## What changed in the conclusions

Three of the reviewer's points changed a conclusion rather than its wording, and we
flag them rather than leaving them for a reader to notice.

1. **Fig. 3f no longer supports a factor of two.** The claim is withdrawn and the
   panel is reported as a demonstration of bias in a cycle-blind fit.
2. **Fig. 8a holidays do beat continuous dosing without a fitness cost**, and the
   reason is a drug effect on the switching rates that was previously stated only in
   the Methods.
3. **The Fig. 4f prediction is confirmed by one experiment and contradicted by
   another**, and which one applies depends on whether the pretreatment biases the
   switching or only speeds it.
