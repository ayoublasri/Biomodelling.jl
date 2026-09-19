"""Python part of Figure 5: GENIE3 (random-forest importances), imputation (MAGIC, kNN-smoothing),
and AUPR/AUROC metrics for all score matrices. Reads/writes paper/output/."""
import glob, os, sys, warnings
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
warnings.filterwarnings("ignore")
OUT = os.path.join(os.path.dirname(__file__), "..", "output")

def genie3(Y, n_trees=200, seed=0):
    n, g = Y.shape
    M = np.zeros((g, g))
    for j in range(g):
        X = np.delete(Y, j, axis=1)
        rf = RandomForestRegressor(n_estimators=n_trees, max_features="sqrt", random_state=seed, n_jobs=4)
        rf.fit(X, Y[:, j])
        imp = rf.feature_importances_
        M[np.arange(g) != j, j] = imp          # regulator (row) -> target (column)
    return M

def knn_smoothing(Y, k=15, n_pcs=10, seed=0):
    """kNN-smoothing (Wagner, Yan and Yanai 2018), one-step variant on Freeman-Tukey transformed PCA space."""
    S = Y.astype(float)
    ls = S.sum(1, keepdims=True); ls[ls == 0] = 1
    N = S / ls * np.median(ls)
    F = np.sqrt(N) + np.sqrt(N + 1)
    Z = PCA(n_components=min(n_pcs, F.shape[1] - 1), random_state=seed).fit_transform(F)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Z)
    idx = nn.kneighbors(Z, return_distance=False)
    return np.array([S[idx[i]].sum(0) for i in range(S.shape[0])])

def metrics(scores, adj, directed):
    g = adj.shape[0]
    if directed:
        mask = ~np.eye(g, dtype=bool)
        truth = (adj != 0)[mask]; s = scores[mask]
    else:
        iu = np.triu_indices(g, 1)
        truth = ((adj != 0) | (adj.T != 0))[iu]; s = np.maximum(scores, scores.T)[iu]
    return average_precision_score(truth, s), roc_auc_score(truth, s), truth.mean()

adj = pd.read_csv(os.path.join(OUT, "fig5b_adjacency.csv")).values
datasets = {os.path.basename(p)[len("fig5b_data_"):-4]: pd.read_csv(p).values for p in glob.glob(os.path.join(OUT, "fig5b_data_*.csv")) if "cycle_regressed" not in p}
# an extra dataset: concentrations with cell-cycle covariates (gene copy number, age) regressed out, using the
# ground-truth metadata the simulator provides — what a perfect cell-cycle correction could achieve
meta = pd.read_csv(os.path.join(OUT, "fig5b_population_metadata.csv"))
Z = np.log1p(datasets["population_concentration"])
X = np.column_stack([np.ones(len(meta)), meta.copies.values, meta.age.values, meta.age.values ** 2])
beta = np.linalg.lstsq(X, Z, rcond=None)[0]
resid = Z - X @ beta
datasets["population_cycle_regressed"] = np.expm1(resid - resid.min(0))       # back to a non-negative scale
pd.DataFrame(datasets["population_cycle_regressed"], columns=[f"gene_{j+1}" for j in range(Z.shape[1])]).to_csv(os.path.join(OUT, "fig5b_data_population_cycle_regressed.csv"), index=False)
for name in ("population_cycle_regressed",):
    Y = datasets[name]; Zc = np.log1p(Y)
    Cm = np.abs(np.corrcoef(Zc.T)); np.fill_diagonal(Cm, 0)
    pd.DataFrame(Cm, columns=[f"gene_{j+1}" for j in range(Cm.shape[1])]).to_csv(os.path.join(OUT, f"fig5b_scores_pearson_{name}.csv"), index=False)
    R = pd.DataFrame(Y).rank().values; Sm = np.abs(np.corrcoef(R.T)); np.fill_diagonal(Sm, 0)
    pd.DataFrame(Sm, columns=[f"gene_{j+1}" for j in range(Sm.shape[1])]).to_csv(os.path.join(OUT, f"fig5b_scores_spearman_{name}.csv"), index=False)
rows = []
# (b) GENIE3 on every dataset
for name, Y in datasets.items():
    M = genie3(np.log1p(Y))
    pd.DataFrame(M, columns=[f"gene_{j+1}" for j in range(M.shape[1])]).to_csv(os.path.join(OUT, f"fig5b_scores_genie3_{name}.csv"), index=False)
# (c) imputation on the sequenced counts, then all inference methods
Yseq = datasets["sequenced_counts"]
imputed = {"knn_smoothing": knn_smoothing(Yseq)}
try:
    import magic
    op = magic.MAGIC(random_state=0, verbose=False, n_jobs=4)
    ls = Yseq.sum(1, keepdims=True); ls[ls == 0] = 1
    imputed["magic"] = op.fit_transform(np.sqrt(Yseq / ls * np.median(ls)))
except Exception as e:
    print("MAGIC skipped:", e)
for name, Y in imputed.items():
    Z = np.log1p(np.maximum(Y, 0))
    C = np.abs(np.corrcoef(Z.T)); np.fill_diagonal(C, 0)
    pd.DataFrame(C, columns=[f"gene_{j+1}" for j in range(C.shape[1])]).to_csv(os.path.join(OUT, f"fig5b_scores_pearson_imputed_{name}.csv"), index=False)
    M = genie3(Z)
    pd.DataFrame(M, columns=[f"gene_{j+1}" for j in range(M.shape[1])]).to_csv(os.path.join(OUT, f"fig5b_scores_genie3_imputed_{name}.csv"), index=False)
# metrics for every score matrix
for p in sorted(glob.glob(os.path.join(OUT, "fig5b_scores_*.csv"))):
    tag = os.path.basename(p)[len("fig5b_scores_"):-4]
    method, dataset = tag.split("_", 1)
    S = pd.read_csv(p).values
    directed = method == "genie3"
    aupr, auroc, base = metrics(S, adj, directed)
    rows.append((method, dataset, directed, aupr, auroc, base))
    print(f"{method:9s} {dataset:32s} AUPR={aupr:.3f} AUROC={auroc:.3f} (random AUPR={base:.3f})")
pd.DataFrame(rows, columns=["method", "dataset", "directed", "aupr", "auroc", "random_aupr"]).to_csv(os.path.join(OUT, "fig5bc_metrics.csv"), index=False)
# (d) perturbation prediction baselines
kd = pd.read_csv(os.path.join(OUT, "fig5d_knockdowns.csv"))
rows = []
for j, grp in kd.groupby("knocked_gene"):
    grp = grp[grp.gene != j]
    true = grp.log2fc_true.values
    lfc_j = kd[(kd.knocked_gene == j) & (kd.gene == j)].log2fc_true.values[0]
    pred_corr = grp.correlation_with_knocked.values * lfc_j
    pred_zero = np.zeros_like(true)
    sst = ((true - true.mean()) ** 2).sum()
    r2 = lambda pred: 1 - ((true - pred) ** 2).sum() / sst if sst > 0 else np.nan
    direct = grp.direct_edge.values != 0
    rows.append((j, lfc_j, r2(pred_zero), r2(pred_corr), np.abs(true[direct]).mean() if direct.any() else np.nan,
                 np.abs(true[~direct]).mean(), direct.sum(), (np.abs(true) > 0.5).sum()))
pd.DataFrame(rows, columns=["knocked_gene", "log2fc_self", "r2_zero_baseline", "r2_correlation_baseline", "mean_abs_lfc_direct_targets",
                            "mean_abs_lfc_indirect", "n_direct_targets", "n_genes_changed"]).to_csv(os.path.join(OUT, "fig5d_summary.csv"), index=False)
print("benchmarks_python done")
