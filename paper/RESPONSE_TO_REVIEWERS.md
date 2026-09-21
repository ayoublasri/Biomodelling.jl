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
Table 5, and the population layer passes.

Along single lineages the mother-machine control makes the recorded cells independent,
so the Kolmogorov-Smirnov test applies as it stands: all four models match the exact law
of their mode, with distances of 0.0023 to 0.0054 against 99% critical values of 0.0081
to 0.0115. A population snapshot is a different matter, and we say so rather than
quoting a test that does not apply: its cells share ancestors, so the independent
critical value understates the true one. The population mode is therefore judged on the
agreement of its replicate means with the exact mean, and over the four models those sit
within 2.5 standard errors; 18 of the 20 individual samples fall below the independent
critical value as well, the two exceptions being replicates of the model with a
deterministic cycle, where cells within a clone are synchronised and the within-sample
correlation is strongest.

Scored against the exact law of the *other* mode, the same samples give distances at
least six times larger, so the comparison has ample power to tell the two settings
apart. The mean molecule number is 13.33 along a lineage against 10.00 in a population
for constitutive production, 13.33 against 10.00 for bursts of four, 8.89 against 6.67
for the two-state promoter, and 25.64 against 23.96 for volume-scaled synthesis with
replication.

The step-size bias is quantified as a theoretical quantity rather than inferred from
noisy samples: the distance between the stationary law of the scheme and that of the
continuous-time model is 0.0004 at a step of 1/512 of the interdivision time and 0.0029
at 1/64, first order in the step and below the critical value at these sample sizes
throughout. In the model with gene replication the cycle-phase distribution the
population develops matches the predicted two-level step to a total variation distance
of 0.020 over twenty bins, so the snapshot weighting emerges from the branching dynamics
rather than being imposed on it.

The Methods describe the protocol, and the closed-form lineage and population means and
variances are now checked on every run of the test suite
(`test/test_exact_population.jl`). The theory itself was written twice, as a matrix
solution of the scheme and as closed-form factorial moments of the continuum limit, and
the two agree to six digits.

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

The matched arm confirms the reviewer's reading and then strengthens the conclusion the
comparison was drawn for. Raising the maximal hazard restores the cycle-blind dose
response almost exactly: the decay rate runs from +0.025 to +0.049 per time unit over
doses 0.5 to 2, against +0.023 to +0.050 with no gate at all and -0.001 to +0.023 at the
unmatched gated hazard. So the shallower response under the gate was indeed the lower
average hazard.

The kin result improves. The question the gate raises is whether it can manufacture the
sister fate correlations without any heritable state, since sisters are born together
and so share a cycle phase. At matched mean hazard it leaves no excess sister
concordance at all in the fast-switching control, which has no usable memory (-0.000,
against +0.008 unmatched), while the memory gene keeps +0.051. The kin signature
therefore reads expression memory rather than the cycle, and that conclusion does not
rest on the gate also killing less. In the melanoma study the best schedule is unchanged
in all three mechanisms at the matched hazard, and under partial protection continuous
dosing still outlasts the intermittent schedule, which is the trial-consistent
ordering.

### 7. Seed replication on the schedule scans, and Supplementary Table 1

Both done. Every point of the Fig. 4h schedule scan and of its cycle-gated counterpart
is now the mean of five independent seeds, and Fig. 4h plots those means. The reviewer's
concern is borne out: the gap between the best schedule and the next best is 0.0003 to
0.0006 per time unit against a standard error of 0.0001 to 0.0008 on the difference, so
the ranking of the top two is resolved in only one of the three mechanisms. The text now
says that, and says that elsewhere the scan identifies a region of good schedules rather
than a single best one. The qualitative conclusions are unaffected: continuous dosing at
the highest dose is still best without a fitness cost, and dose 1 still beats dose 2 for
drug-induced tolerance.

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
added to the Discussion, with the trial's primary numbers stated: the progression-free
survival advantage of continuous dosing (median 9.0 against 5.5 months, P = 0.063
against the pre-specified two-sided alpha of 0.2) did not carry through to overall
survival, a secondary end point the trial was not powered for and on which the two arms
were equal at a median of 29.2 months; and the resistant xenograft tumours were drug-addicted, regressing on withdrawal
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

---

# Response to the third-round report

Every point is addressed below. Two of the referee's findings could not be reproduced
against the sources, and we set out the evidence rather than make a change that would
introduce an error; three others we carried further than asked, because the underlying
defect turned out to be larger than the symptom reported.

## Blocking

**B1. The stale cycle-gate paragraph.** Confirmed, and the cause was not the one either
possibility in the report anticipated. The gated scan *was* re-run on the same five
seeds; what was stale was the cycle-blind *reference* it was compared against.
`fill_numbers.py` loaded `fig4f_schedules.csv`, the older single-seed table, while the
gated arm read the five-seed rerun, so the apparent release-period difference was a seed
artefact. Against the matching five-seed scan no release period moves at all.

The generator now loads the five-seed table and reports the real comparison: continuous
dosing wins in all three mechanisms with the gate as without it, at the same dose in
three of three. The gate raises the long-term growth rate at all sixty points of the
scan (+0.0026 to +0.0155 per time unit) and perturbs the order below the winner
(Spearman correlation 0.977 to 0.998), but changes neither which schedule wins nor the
intermediate-dose optimum under drug-induced tolerance. We also state that the winner's
margin over the best holiday exceeds twice its standard error in five of the six
comparisons, the exception being pre-existing tolerance without the gate
(+0.0004 against a standard error of 0.0005), so the claim is hedged to the same degree
as the neighbouring seed note. Supplementary Note 13 and the Results are generated from
the same text, so both are fixed at once; Supplementary Fig. 5's caption covers only the
decay, concordance and melanoma panels and never carried the claim.

**B2. The pre-resistant fraction.** Confirmed; option (b) adopted. The Results now say
the melanoma study does *not* sit in the traced regime, give the value (1:200, the
Shaffer per-marker frequency), state that it is 5- to 50-fold higher, and note that a
draw of the study's 2,000 founders would contain no pre-resistant cell at all with
probability 0.8 at 1 in 10,000. The direction of the bias is named: among the arms that
do lose control, control is lost earlier than the traced frequency would give. We
restricted that statement to those arms deliberately, because sixteen of the thirty-six
arms are censored at 52 weeks and their reported values would be unchanged at a lower
fraction. We do not claim statistical tractability as the *reason* the value was chosen,
since nothing in the record establishes a motive.

**B3. The temozolomide half-life.** Confirmed, and the correction needed to go further
than relabelling. Both 2.1 h and 1.8 h are plasma elimination half-lives, so the
contrastive framing itself was wrong, not just the "cerebrospinal-fluid" label: replacing
it with "systemic" would have preserved the false implication that 2.1 h is not the
plasma value, and would have contradicted the Methods sentence the report ruled correct.
All three locations now describe 2.1 h as a plasma value from Ostermann's
three-compartment population model of plasma and cerebrospinal fluid, and 1.8 h as an
independent single-dose estimate. The genuinely CSF-specific quantity, AUC_CSF equal to
20% of AUC_plasma, is added to Supplementary Table 4 as its own row. The Methods
parenthesis attributed that same 20% to a brain-tumour compartment; Ostermann measured
cerebrospinal fluid, not tumour tissue, so that sentence is corrected too.

**B15. Affiliation.** Outstanding, and deliberately so. The mechanism is settled: pandoc's
default LaTeX template carries no affiliation field, so it has to ride on the title-page
line beneath the author, which we verified renders correctly in both PDF and DOCX. What
we have not done is assert an affiliation on the author's behalf. "Independent
researcher" is a positive claim about institutional status that nothing in the record
establishes, and an ORCID cannot be invented at all. Both are supplied by the author
before submission.

## Not reproduced

**Checklist 11. The Piho & Thomas attribution.** The report has the two papers the wrong
way round. The sentence quoted as the title of the PLoS Computational Biology paper — "a
finite state projection approach to analyse gene expression and division distributions
and infer selection from single-cell data in mother machines and lineage trees" — is a
sentence from the *abstract* of the Science Advances 2024 paper. The PLoS Computational
Biology paper is titled "AgentBasedModeling.jl: a tool for stochastic simulation of
structured population dynamics" and contains no finite-state-projection inference claim.
Table 1's caption already attaches the inference claim to the Science Advances paper and
the tool itself to the PLoS paper, which is correct as it stands. Both bibliography
entries match their records exactly, so nothing is changed.

**B5, the sqrt(5) claim.** Not reproduced. The code already divides the spread between
the five replicate means by the square root of five: `fill_numbers.py` computes
`g.mean_sim.std(ddof=1) / np.sqrt(len(g))` for the population rows and uses the
within-sample error only for the single-sample lineage rows. Recomputing from
`fig9_validation.csv` gives deviations of -0.71, -0.70, -2.45 and -0.08 standard errors,
which are the values the table prints. The two-state promoter is 2.45 standard errors
low, not 5.5; on Student's t with four degrees of freedom that is p = 0.07.

What misled the report is real, and was our error: the Supplementary Table 5 caption
described the quoted error as "the spread between them" rather than as the standard error
of the average. Read that way, the sqrt(5) objection follows exactly. The caption now
states the statistic explicitly.

The report is also understating one thing, which we have taken up. It is not three of
four population rows that are low but **four of four**, by 0.03% to 0.70%. No single
deviation is resolved, but the common sign deserved an account, so we tested the three
candidate causes and excluded all of them: the discretisation cannot be responsible
because each population sample is scored against the stationary law of the *scheme* at
the simulated step, and solving that law from 1/64 to 1/512 of the interdivision time
moves the exact population mean by less than one part in 10^9; the subsampling cap draws
uniformly without replacement from a snapshot whose division and thinning are independent
of the molecule number, so the retained cells are exchangeable; and twelve generations of
burn-in leave a transient far below the observed offset. The text now says that an offset
of this size is neither established nor excluded, rather than asserting agreement.

## The numbered checklist

**B6.** Corrected to "a hundred time units, five cell-cycle times".
**B7.** The colour scale is Fig. 7e. The panel letter is fixed and the Fig. 7e legend,
which described the scale as a distance to the training data, now says it is the
objective the optimiser minimised, which adds the memory prior and is not the
seed-averaged error quoted in the text.
**B8.** Corrected, and not in the direction proposed. Reporting the fast control's
enrichment as "small (1.8- and 2.5-fold)" would assert an enrichment the data do not
support: those folds rest on one sister pair and two cousin pairs in which both lineages
survived, and the control's excess concordance is +0.000. It would also have contradicted
the cycle-gate paragraph six lines later, which already reports that +0.000. The text now
says what the folds rest on and that the concordance sits at the independent expectation.
**B9.** Corrected, and narrowed. Only the *training* pair of the killing-mechanism
comparison is single-seed; its held-out pair is averaged over three seeds. Both are
quoted in the main text, so the exception is stated there rather than sited in the
supplement.
**B10.** Both references to an earlier version are removed.
**B11.** The Fig. 3f values are now identified as scaled, with the unscaled rate given
and its missing time dimension restored.
**B12.** Corrected. The interval was a normal approximation; on three seeds it is now
Student's t on two degrees of freedom, and the method is stated. The interval widens from
+0.000 to +0.043 to -0.026 to +0.069, so it contains zero — which supports rather than
weakens the claim that third-cousin correlations are absent.
**B13.** Reconciled, and computed rather than asserted: the scan gives 2.3 weeks over
three seeds, the main run 2.7 weeks over four with a standard error of 0.2 weeks, and
differencing the rounded 23 and 26 gives 3. Supplementary Note 15 now prints all three
and says they are one result.
**B14.** The Introduction now separates the statistical and reference-based simulators
(Splatter, scDesign3, GRouNdGAN) from the network-based ones (SymSim, SERGIO, dyngen,
scMultiSim), with each citation attached to the right tool. We did not restrict
trajectory generation to the network-based group, because Splatter, scDesign3 and SymSim
all generate trajectories. The same blanket claim survived in the **abstract**, which the
report did not flag and which contradicted the Introduction and Discussion after the
second round; it is corrected, and the abstract is 199 words with no citations.

## Bibliography

Ahlmann-Eltze now carries 22(8), 1657-1661; scMultiSim 22(5), 982-993; and the
DifferentialEquations.jl entry was typed `@misc`, which pandoc maps to an empty CSL type
so that `nature.csl` never reached its journal branch — it is now `@article` and renders
as J. Open Res. Softw. 5, 15 (2017).

The duplicated DOIs were a style bug, not a data one. `nature.csl` prints a doi.org URL
for any entry without a volume and then prints the DOI again through its `access` macro.
Stripping the `doi` fields would have hidden the symptom at the cost of the machine-readable
metadata, so the macro is fixed instead: it now suppresses the second form whenever a DOI
is present. All seven affected entries render once and keep their `doi` fields.

The Corigliano article number could not be confirmed. Every primary source is blocked by
this environment's egress policy (doi.org, crossref, OpenAlex and journals.aps.org all
return 403 on CONNECT), and search snippets are not an acceptable source for a citation.
Following the report's own fallback, the entry cites the DOI and prints no volume or
article number; a check from an unrestricted network would settle it.

## SWOG S1320

Corrected, and the claim withdrawn rather than re-attributed. Overall survival was a
median of 29.2 months in both arms, so "favoured the intermittent arm" was not supportable
in the form it was written. The Discussion now gives the primary result with the trial's
own threshold, which matters: S1320 pre-specified a two-sided alpha of 0.2 and reported
80% confidence intervals, so P = 0.063 *is* the significant progression-free survival
advantage the sentence describes, and quoting it without the threshold would have made
the manuscript look wrong. Overall survival is labelled as the secondary end point the
trial was not powered for. We did not substitute the prior-immunotherapy subgroup result,
because its numbers could not be verified from this environment; an author with access
may add it, flagged as post hoc.

## Repository

Every file named in the report exists: `paper/run_all.jl`, `paper/scripts/fig7_calibration.jl`,
`paper/scripts/fig9_validation.jl`, `paper/scripts/exact_solutions.py`,
`paper/manuscript/fill_numbers.py`, `test/test_exact_population.jl`,
`docs/src/migration.md` and a genuine MIT `LICENSE`. GitHub Actions CI is present
(`CI.yml`, `Docs.yml`, `CompatHelper.yml`, `TagBot.yml`, `paper-compute.yml`) and
`Project.toml` declares version 2.0.0. Two real defects in the reproduction recipe were
found and fixed: `run_all.jl` does not regenerate Supplementary Fig. 6 or the cell-cycle
panels, so the opt-in invocations are now named in Methods; and Methods claimed the
four-thread *timings* are in `paper/README.md` when what the README gives is the
four-thread *invocation*.

Two things remain for the author and cannot be done here: **there is no v2.0.0 tag or
release** (`git tag -l` is empty), and **no Zenodo DOI** appears in the Data or Code
availability statements. The Data availability sentence asserts that the generated tables
are archived with the tagged v2.0.0 release, which becomes true only once that release is
cut. We did not insert a placeholder DOI, because it would print as a broken identifier
in the submitted PDF.
