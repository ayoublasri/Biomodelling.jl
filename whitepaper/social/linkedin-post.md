# LinkedIn post

Image: `linkedin-card.png` (1200 x 1200, rendered at 2x). Rebuild it with
`python build_card.py` and a headless Chromium screenshot of `linkedin-card.html`.

---

Why does the same drug schedule win in one study and lose in another?

In the SWOG S1320 melanoma trial, continuous BRAF/MEK inhibition beat intermittent dosing: 9.0 against 5.5 months of median progression-free survival. In patient-derived xenografts, intermittent dosing was the schedule that worked. Both results are real.

I rebuilt my simulation framework to find out what separates them.

Biomodelling.jl 2.0 simulates single cells that grow, divide, inherit their gene-expression states and respond to drugs, one stochastic reaction at a time. Drug tolerance is not assumed. It emerges from a single gene whose promoter switches more slowly than the cell divides, so a rare protected state is passed to daughters.

In the model, one property of the resistant cells decides the answer:

• Resistance is free. Holidays only let sensitive cells regrow and be killed again. Continuous dosing wins.
• Resistance costs growth. Holidays let sensitive cells outcompete the resistant ones. Adaptive therapy keeps control for a full year on roughly a third of the cumulative dose.
• Resistance is partial. Holidays release the resistant cells. Continuous dosing wins again, which is the ranking the trial reported.

That property is measurable in a dish, before a schedule is chosen.

The model earns the right to say this. It was calibrated on published time-lapse data of cisplatin-treated cells, then tested on a concentration it never saw, on the correlation of fates between sisters and cousins, and on the observation that single-cell death times barely shift with dose while population kill rates change several fold. A built-in optimiser then searches schedule space instead of guessing.

Everything is open. Code under MIT, the white paper, and every figure reproducible from fixed seeds:
github.com/ayoublasri/Biomodelling.jl

One more thing worth saying plainly: version 2.0 and the white paper were authored with Claude (Anthropic) working under my direction. The rewrite, the simulations, the calibration and the figures came out of that collaboration, and I verified the results.

What would you want to measure first, if the schedule depended on it?

#SingleCell #DrugDiscovery #ComputationalBiology #CancerResearch #Julia #OpenScience #SystemsBiology
