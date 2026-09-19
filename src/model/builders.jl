# ---------------------------------------------------------------------------
# Convenience model builders
# ---------------------------------------------------------------------------

"""
    birth_death_model(; k, γ, species=:X, volume_scaled=false)

Constitutive production `∅ → X` (rate `k`, ∝ V if `volume_scaled`) and first-order
degradation `X → ∅` (rate `γ`). Stationary distribution: Poisson(k/γ) at V = 1.
"""
function birth_death_model(; k, γ, species::Symbol=:X, volume_scaled::Bool=false)
    ReactionModel([
        Reaction("production", [], [species], MassAction(:k); volume = volume_scaled ? :proportional : :none),
        Reaction("degradation", [species], [], MassAction(:γ)),
    ]; params = (k = k, γ = γ))
end

"""
    telegraph_model(; k_on, k_off, k_tx, k_dm, k_tl=nothing, k_dp=nothing,
                    gene=:G, mrna=:mRNA, protein=:protein, volume_scaled=true)

Two-state promoter (`G_off ⇄ G_on`), transcription from the active state and
mRNA degradation; optionally translation and protein degradation. The promoter
species form a promoter group. Stationary mRNA distribution at fixed volume:
Beta-Poisson (Peccoud and Ycart 1995), see [`telegraph_pmf`](@ref).
"""
function telegraph_model(; k_on, k_off, k_tx, k_dm, k_tl=nothing, k_dp=nothing,
                         gene::Symbol=:G, mrna::Symbol=:mRNA, protein::Symbol=:protein,
                         volume_scaled::Bool=true)
    Goff = Symbol(gene, "_off"); Gon = Symbol(gene, "_on")
    rx = Reaction[
        Reaction("activation", [Goff], [Gon], MassAction(:k_on)),
        Reaction("inactivation", [Gon], [Goff], MassAction(:k_off)),
        Reaction("transcription", [Gon], [Gon, mrna], MassAction(:k_tx); volume = volume_scaled ? :proportional : :none),
        Reaction("mRNA degradation", [mrna], [], MassAction(:k_dm)),
    ]
    params = Dict{Symbol,Float64}(:k_on => k_on, :k_off => k_off, :k_tx => k_tx, :k_dm => k_dm)
    if k_tl !== nothing
        k_dp === nothing && throw(ArgumentError("k_dp is required with k_tl"))
        push!(rx, Reaction("translation", [mrna], [mrna, protein], MassAction(:k_tl); volume = :none))
        push!(rx, Reaction("protein degradation", [protein], [], MassAction(:k_dp)))
        params[:k_tl] = k_tl; params[:k_dp] = k_dp
    end
    ReactionModel(rx; params = params, promoters = [[Goff, Gon]])
end

"""
    two_stage_model(; k_tx, k_dm, k_tl, k_dp, mrna=:mRNA, protein=:protein, volume_scaled=false)

Constitutive transcription, translation and degradation (the classic two-stage
model). With short-lived mRNA the protein distribution approaches a negative
binomial (Shahrezaei and Swain 2008).
"""
function two_stage_model(; k_tx, k_dm, k_tl, k_dp, mrna::Symbol=:mRNA, protein::Symbol=:protein, volume_scaled::Bool=false)
    ReactionModel([
        Reaction("transcription", [], [mrna], MassAction(:k_tx); volume = volume_scaled ? :proportional : :none, copy_number = true),
        Reaction("mRNA degradation", [mrna], [], MassAction(:k_dm)),
        Reaction("translation", [mrna], [mrna, protein], MassAction(:k_tl); volume = :none),
        Reaction("protein degradation", [protein], [], MassAction(:k_dp)),
    ]; params = (k_tx = k_tx, k_dm = k_dm, k_tl = k_tl, k_dp = k_dp))
end

"""
    bursty_protein_model(; a, b, γ, protein=:protein)

Protein produced in geometric bursts of mean size `b` arriving at rate `a`, with
degradation `γ`, implemented as a `Custom` kinetics burst channel. The stationary
distribution is negative binomial with shape `a/γ` and success probability
`1/(1+b)` (Friedman et al. 2006, discrete version).
"""
function bursty_protein_model(; a, b, γ, protein::Symbol=:protein)
    # Bursts are emulated by an mRNA that is translated many times before decay:
    # mRNA lifetime → 0 with fixed mean burst size b = k_tl / k_dm.
    k_dm = 200.0 * γ
    ReactionModel([
        Reaction("burst", [], [:_m], MassAction(:a); volume = :none),
        Reaction("burst end", [:_m], [], MassAction(:k_dm)),
        Reaction("translation", [:_m], [:_m, protein], MassAction(:k_tl); volume = :none),
        Reaction("degradation", [protein], [], MassAction(:γ)),
    ]; params = (a = a, k_dm = k_dm, k_tl = b * k_dm, γ = γ))
end
