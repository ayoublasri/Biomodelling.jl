# ---------------------------------------------------------------------------
# Population control
# ---------------------------------------------------------------------------
abstract type PopulationControl end

"""
    ConstantN()

Keep the number of cells constant: each second daughter replaces a uniformly
random cell of the population (Moran-like replacement, the behaviour of the
original `exponential_growth`). Cells killed by drug are removed without
replacement.
"""
struct ConstantN <: PopulationControl end

"""
    FreeGrowth(; max_cells=100_000)

Branching population: every daughter is kept. When the population exceeds
`max_cells` it is subsampled uniformly and the true size is tracked through a
scaling factor (see `popsize`).
"""
struct FreeGrowth <: PopulationControl
    max_cells::Int
end
FreeGrowth(; max_cells::Integer=100_000) = FreeGrowth(Int(max_cells))

"""
    LogisticGrowth(K; max_cells=100_000)

Free growth with an additional death hazard `λ N / K` (carrying capacity `K`).
"""
struct LogisticGrowth <: PopulationControl
    K::Float64
    max_cells::Int
end
LogisticGrowth(K::Real; max_cells::Integer=100_000) = LogisticGrowth(Float64(K), Int(max_cells))

# ---------------------------------------------------------------------------
# Settings and result
# ---------------------------------------------------------------------------

"""
    PopulationSettings(; dt=0.1, kernel=nothing, growth=ExponentialGrowth(0.03),
                        size_control=Sizer(2.0), partitioning=BinomialPartition(),
                        replication=nothing, control=ConstantN(), background_death=0.0,
                        record_every=1, track_lineage=true, randomize_promoters=true,
                        threads=true)

Options of [`simulate_population`](@ref). `dt` is the global update and recording
step (volume is held constant within a step); `kernel` defaults to
`HybridSSATau(dt)`.
"""
Base.@kwdef struct PopulationSettings
    dt::Float64 = 0.1
    kernel::Union{Nothing,AbstractKernel} = nothing
    growth::GrowthModel = ExponentialGrowth(0.03)
    size_control::SizeControl = Sizer(2.0)
    partitioning::Partitioning = BinomialPartition()
    replication::Union{Nothing,Replication} = nothing
    control::PopulationControl = ConstantN()
    background_death::Float64 = 0.0
    record_every::Int = 1
    track_lineage::Bool = true
    randomize_promoters::Bool = true
    threads::Bool = true
end

"""
    PopulationResult

Output of [`simulate_population`](@ref). Per recorded time `t[i]`: `counts[i]`
(cells × species), `volume[i]`, `ids[i]`, `age[i]`, `generation[i]`, `clone[i]`,
`copies[i]`, `perturbed[i]`; `popsize[i]` is the (scaled) number of cells and
`dose[i]` the dose applied during the step. `lineage` is the [`LineageTable`](@ref).
"""
struct PopulationResult
    t::Vector{Float64}
    counts::Vector{Matrix{Int}}
    volume::Vector{Vector{Float64}}
    ids::Vector{Vector{Int}}
    age::Vector{Vector{Float64}}
    generation::Vector{Vector{Int}}
    clone::Vector{Vector{Int}}
    copies::Vector{Vector{Int}}
    perturbed::Vector{Vector{Bool}}
    popsize::Vector{Float64}
    dose::Vector{Float64}
    lineage::LineageTable
    model::ReactionModel
    settings::PopulationSettings
    p::Vector{Float64}
end

Base.show(io::IO, r::PopulationResult) = print(io, "PopulationResult(", length(r.t), " records, final N=",
    round(Int, r.popsize[end]), ", ", nspecies(r.model), " species, ", r.lineage, ")")

"""
    snapshot(result, i) -> NamedTuple
    snapshot(result; t) -> NamedTuple

Recorded population state at record `i` (or at the record closest to time `t`).
"""
function snapshot(r::PopulationResult, i::Integer)
    (t = r.t[i], counts = r.counts[i], volume = r.volume[i], ids = r.ids[i], age = r.age[i],
     generation = r.generation[i], clone = r.clone[i], copies = r.copies[i], perturbed = r.perturbed[i],
     popsize = r.popsize[i], dose = r.dose[i], species = speciesnames(r.model))
end
snapshot(r::PopulationResult; t::Real) = snapshot(r, argmin(abs.(r.t .- t)))
final_snapshot(r::PopulationResult) = snapshot(r, length(r.t))
popsize(r::PopulationResult) = r.popsize

"""
    concentrations(result, i) -> Matrix{Float64}

Counts divided by cell volume at record `i`.
"""
concentrations(r::PopulationResult, i::Integer) = r.counts[i] ./ r.volume[i]

speciesindex(r::PopulationResult, s) = speciesindex(r.model, s)

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

function _initial_volume(sc::SizeControl, Vb::Float64, target::Float64, rng::AbstractRNG)
    hi = sc isa AgeTimer ? 2Vb : max(target, Vb * 1.01)
    Vb + rand(rng) * (hi - Vb)
end

function _daughter(c::Cell, x::Vector{Int}, V::Float64, t::Float64, next_id::Base.RefValue{Int}, settings::PopulationSettings)
    crng = Xoshiro(rand(c.rng, UInt64))
    target = sample_target(settings.size_control, V, crng)
    λ = sample_growth_rate(settings.growth, crng)
    Cell(next_id[] += 1, c.id, c.clone, c.generation + 1, t, V, V, target, 1, false, λ, x, copy(c.p), crng, c.perturbed, true)
end

function divide!(c::Cell, m::ReactionModel, settings::PopulationSettings, t::Float64, lt::LineageTable,
                 next_id::Base.RefValue{Int}, base_copies::Vector{Int})
    f = partition_fraction(settings.partitioning, c.rng)
    xa = similar(c.x); xb = similar(c.x)
    partition!(xa, xb, c.x, f, m, settings.partitioning, base_copies, c.rng)
    settings.track_lineage && record_end!(lt, c.id, t, :divided, c.x, c.V)
    a = _daughter(c, xa, c.V * f, t, next_id, settings)
    b = _daughter(c, xb, c.V * (1 - f), t, next_id, settings)
    if settings.track_lineage
        record_birth!(lt, a); record_birth!(lt, b)
    end
    a, b
end

function _record!(ts, counts, volume, ids, ages, gens, clones, copies, perturbed, popsizes, doses,
                  t::Float64, cells::Vector{Cell}, N::Float64, d::Float64, ns::Int)
    n = length(cells)
    C = Matrix{Int}(undef, n, ns)
    V = Vector{Float64}(undef, n); I = Vector{Int}(undef, n); A = Vector{Float64}(undef, n)
    G = Vector{Int}(undef, n); K = Vector{Int}(undef, n); Cp = Vector{Int}(undef, n); P = Vector{Bool}(undef, n)
    for (i, c) in enumerate(cells)
        C[i, :] .= c.x
        V[i] = c.V; I[i] = c.id; A[i] = t - c.birth_time; G[i] = c.generation; K[i] = c.clone; Cp[i] = c.copies; P[i] = c.perturbed
    end
    push!(ts, t); push!(counts, C); push!(volume, V); push!(ids, I); push!(ages, A); push!(gens, G)
    push!(clones, K); push!(copies, Cp); push!(perturbed, P); push!(popsizes, N); push!(doses, d)
    nothing
end

function _step_cell!(c::Cell, m::ReactionModel, pp::Vector{Float64}, cp, d::Float64, t::Float64, tn::Float64,
                     dt::Float64, kernel::AbstractKernel, ws::Workspace, settings::PopulationSettings, extra_h::Float64)
    gmult, h = apply_effects!(c, cp, pp, d, t)
    simulate!(c.x, m, c.p, c.V, c.copies, t, tn, kernel, ws, c.rng)
    c.V = grow(settings.growth, c.V, c.growth_rate * gmult, dt)
    if settings.replication !== nothing && !c.replicated &&
       cycle_progress(settings.size_control, c, tn) >= settings.replication.fraction
        replicate!(c, m)
    end
    htot = h + settings.background_death + extra_h
    if htot > 0.0 && rand(c.rng) < -expm1(-htot * dt)
        c.alive = false
    end
    nothing
end

"""
    simulate_population(model, x0, N0, tspan; settings=PopulationSettings(), p=model.p0,
                        perturbation=nothing, rng=Random.default_rng(), V0=nothing) -> PopulationResult

Simulate `N0` initial cells with state `x0` from `tspan[1]` to `tspan[2]`. Each
cell runs the stochastic kinetics of `model` at its current volume, grows,
replicates its genes (optional), divides with partitioning of molecules, and may
die under a [`Perturbation`](@ref). Every cell carries its own random number
generator seeded from `rng`, so results are reproducible regardless of the
number of threads.
"""
function simulate_population(m::ReactionModel, x0::AbstractVector{<:Integer}, N0::Integer, tspan::Tuple{<:Real,<:Real};
                             settings::PopulationSettings=PopulationSettings(), p::AbstractVector{<:Real}=m.p0,
                             perturbation::Union{Nothing,Perturbation}=nothing, rng::AbstractRNG=Random.default_rng(),
                             V0=nothing)
    t0, t1 = float(tspan[1]), float(tspan[2])
    dt = settings.dt
    nsteps = max(1, round(Int, (t1 - t0) / dt))
    kernel = settings.kernel === nothing ? HybridSSATau(dt) : settings.kernel
    pp = Vector{Float64}(p)
    x0v = Vector{Int}(x0)
    length(x0v) == nspecies(m) || throw(ArgumentError("x0 must have $(nspecies(m)) entries"))
    cp = perturbation === nothing ? nothing : compile(perturbation, m)
    base_copies = [sum(x0v[s] for s in g) for g in m.promoter_groups]
    any(==(0), base_copies) && throw(ArgumentError("every promoter group needs at least one copy in x0"))
    lt = LineageTable()
    next_id = Ref(0)
    sc = settings.size_control
    Vb = typical_birth_volume(sc)
    cells = Vector{Cell}(undef, N0)
    for i in 1:N0
        crng = Xoshiro(rand(rng, UInt64))
        x = copy(x0v)
        settings.randomize_promoters && randomize_promoters!(x, m, base_copies, crng)
        target = sample_target(sc, Vb, crng)
        V = V0 === nothing ? _initial_volume(sc, Vb, target, crng) : Float64(V0[i])
        λ = sample_growth_rate(settings.growth, crng)
        birth = sc isa AgeTimer ? t0 - rand(crng) * target : t0
        c = Cell(next_id[] += 1, 0, i, 0, birth, V, Vb, target, 1, false, λ, x, copy(pp), crng, false, true)
        cells[i] = c
        settings.track_lineage && record_birth!(lt, c)
    end
    ns = nspecies(m)
    ts = Float64[]; counts = Matrix{Int}[]; volume = Vector{Float64}[]; ids = Vector{Int}[]
    ages = Vector{Float64}[]; gens = Vector{Int}[]; clones = Vector{Int}[]; copiesv = Vector{Int}[]
    perturbed = Vector{Bool}[]; popsizes = Float64[]; doses = Float64[]
    scale = 1.0
    nth = settings.threads ? Threads.nthreads() : 1
    wss = [Workspace(m) for _ in 1:nth]
    genes_assigned = cp === nothing || isempty(cp.genes)
    _record!(ts, counts, volume, ids, ages, gens, clones, copiesv, perturbed, popsizes, doses,
             t0, cells, N0 * scale, cp === nothing ? 0.0 : dose(cp.schedule, t0), ns)
    for step in 1:nsteps
        t = t0 + (step - 1) * dt
        tn = t + dt
        d = cp === nothing ? 0.0 : dose(cp.schedule, t)
        if !genes_assigned && t >= minimum(g.t_start for g in cp.genes)
            for c in cells
                c.perturbed = rand(c.rng) < cp.genes[1].fraction
            end
            genes_assigned = true
        end
        N = length(cells)
        extra_h = settings.control isa LogisticGrowth ? settings.growth.rate * (N * scale) / settings.control.K : 0.0
        if settings.threads && nth > 1 && N > 1
            let cells = cells
                @sync for (k, chunk) in enumerate(_chunks(N, nth))
                    Threads.@spawn begin
                        local tws = wss[k]
                        for i in chunk
                            local tc = cells[i]
                            tc.alive && _step_cell!(tc, m, pp, cp, d, t, tn, dt, kernel, tws, settings, extra_h)
                        end
                    end
                end
            end
        else
            for i in 1:N
                c = cells[i]
                c.alive && _step_cell!(c, m, pp, cp, d, t, tn, dt, kernel, wss[1], settings, extra_h)
            end
        end
        # deaths
        if any(c -> !c.alive, cells)
            alive = Cell[]
            for c in cells
                if c.alive
                    push!(alive, c)
                else
                    settings.track_lineage && record_end!(lt, c.id, tn, :died, c.x, c.V)
                end
            end
            cells = alive
        end
        # divisions
        newborns = Cell[]
        for i in eachindex(cells)
            c = cells[i]
            if should_divide(sc, c, tn)
                a, b = divide!(c, m, settings, tn, lt, next_id, base_copies)
                cells[i] = a
                push!(newborns, b)
            end
        end
        # population control
        if settings.control isa ConstantN
            for b in newborns
                isempty(cells) && break
                j = rand(rng, 1:length(cells))
                settings.track_lineage && record_end!(lt, cells[j].id, tn, :removed, cells[j].x, cells[j].V)
                cells[j] = b
            end
        else
            append!(cells, newborns)
            maxc = settings.control.max_cells
            if length(cells) > maxc
                keep = sort!(sample(rng, 1:length(cells), maxc; replace=false))
                scale *= length(cells) / maxc
                if settings.track_lineage
                    keepset = Set(keep)
                    for (i, c) in enumerate(cells)
                        i in keepset || record_end!(lt, c.id, tn, :removed, c.x, c.V)
                    end
                end
                cells = cells[keep]
            end
        end
        if step % settings.record_every == 0 || step == nsteps
            _record!(ts, counts, volume, ids, ages, gens, clones, copiesv, perturbed, popsizes, doses,
                     tn, cells, length(cells) * scale, d, ns)
        end
        isempty(cells) && break
    end
    PopulationResult(ts, counts, volume, ids, ages, gens, clones, copiesv, perturbed, popsizes, doses, lt, m, settings, pp)
end
