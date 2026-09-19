# Synthetic single-cell data

True molecule counts of a snapshot are turned into measurements with an
observation model.

```julia
sn = sample_cells(res; n = 2000, rng = Xoshiro(2))          # 2000 random cells of the last record
seq = sequence(sn.counts, SeqProtocol(capture = 0.15, capture_cv = 0.3, depth = 5000.0, dropout = 0.02))
seq.Y            # observed counts (cells × genes)
seq.capture      # per-cell capture efficiency
```

Because counts scale with cell volume, library size is correlated with the
true volume of each cell (`sn.volume`), which reproduces the transcriptome-size
confound of real data.

Other modalities:

```julia
smfish(sn.counts; efficiency = 0.95)
timelapse(res, :protein; interval = 1.0, σ = -0.1, n_lineages = 20)   # reporter traces along lineages
```

Export:

```julia
write_counts_csv("counts.csv", seq.Y, sn.species, sn.ids)
write_metadata_csv("meta.csv", sn)
write_lineage_csv("lineage.csv", res.lineage)
using HDF5
write_h5ad("data.h5ad", seq.Y; genes = sn.species, cells = sn.ids,
           obs = (volume = sn.volume, age = sn.age, generation = sn.generation, clone = sn.clone),
           obsm = (X_true = Float64.(sn.counts),))
```

The `.h5ad` file opens with `anndata.read_h5ad` in Python.

Random networks with known ground truth:

```julia
model, adjacency = random_grn(100; n_activations = 120, n_inhibitions = 60, telegraph = true, topology = :scale_free)
```
