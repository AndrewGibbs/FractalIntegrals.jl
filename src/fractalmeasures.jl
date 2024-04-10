abstract type Fractal end

abstract type AbstractAttractor end

struct HomogenousAttractor{ T<:Real,
                            S<:AbstractArray{<:AbstractSimilarity},
                            G<:AbstractArray{<:AbstractInvariantMap}
                            } <: AbstractAttractor
    ifs::S
    diam::T
    d::T
    n::Int64
    connectedness::Matrix{Bool}
    symmetries::G
    ρ::T
end

struct Attractor{   T<:Real,
                    S<:AbstractArray{<:AbstractSimilarity},
                    G<:AbstractArray{<:AbstractInvariantMap},
                    } <: AbstractAttractor
    ifs::S
    diam::T
    d::T
    n::Int64
    connectedness::Matrix{Bool}
    symmetries::G
end

function ifs_map!(  Sx::AbstractVector{<:T},
                    S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{<:T}) where T
    # M = length(S)
    N = length(x)
    # y = zeros(T, M*N)
    for m in eachindex(S)
        Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
    end
end

function ifs_map(   S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{T}) where T
    N = length(x)
    Sx = zeros(T, length(S)*N)
    for m in eachindex(S)
        @inbounds Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
    end
    return Sx
end

(Γ::AbstractAttractor)(x::AbstractVector) = ifs_map(Γ.ifs, x)

function get_subattractor_elements(Γ::AbstractAttractor, index::AbstractVector{<:Integer})
    # get new measure and diameter. First initialise:
    new_diam = Γ.diam
    new_ifs = Γ.ifs
    new_symmetries = Γ.symmetries

    for m = index[end:-1:1]
        new_diam *= Γ.ifs[m].ρ
        new_ifs = simcompifs(Γ.ifs[m], new_ifs)
        new_symmetries = simcompsymmetries(Γ.ifs[m], new_symmetries)
    end

    return new_diam, new_ifs, new_symmetries
end

function get_subattractor(Γ::HomogenousAttractor, m::AbstractVector{<:Integer})
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, m)
    return HomogenousAttractor(new_ifs,
            new_diam,
            Γ.d,
            Γ.n,
            Γ.connectedness,
            new_symmetries,
            Γ.ρ)
end

function get_subattractor(Γ::A, m::AbstractVector{<:Integer}) where A<:Attractor
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, m)
    return A(new_ifs,
            new_diam,
            Γ.d,
            Γ.n,
            Γ.connectedness,
            new_symmetries)
end

abstract type AbstractInvariantMeasure{T<:Real, B, V<:AbstractVector{T}, A<:AbstractAttractor} end

# I'm wondering if its actuall worth defining different types of measures.
# But I think that HausdorffMeasure{HomogenousAttractor} will be quite useful
# also for symmetries, e.g. Koch is inhomogenous but has symmetries

struct HausdorffMeasure{T <: Real,
                        B, # probably want this to match the translation in A
                        V <: AbstractVector{T},
                        A <: AbstractAttractor,
                        G <: AbstractArray{<:AbstractInvariantMap}
                        } <: AbstractInvariantMeasure{T, B, V, A}
    supp::A
    barycentre::B
    suppmeasure::T
    weights::V
    symmetries::G
end

# make natural conversions of attractors to Hausdorff measures
# HausdorffMeasure(Γ::AbstractAttractor) = HausdorffMeasure(Γ, ...)

# note the situation with symmetries - they need not be the same for an attractor and a measure.
# But for HausdorffMeasure, the symmetries should be automatically inhereted from the attractor.

struct InvariantMeasure{T <: Real,
                        B, # probably want this to match the translation in A
                        V <: AbstractVector{T},
                        A <: AbstractAttractor,
                        G <: AbstractArray{<:AbstractInvariantMap}
                        } <: AbstractInvariantMeasure{T, B, V, A}
    supp::A
    barycentre::B
    suppmeasure::T
    weights::V
    symmetries::G
end

function get_barycentre(sims::AbstractVector{<:AbstractSimilarity},
                        weights::AbstractVector{<:Real})
    M = length(sims)
    divisor = IdMat
    vec_sum = zero(sims[1].δ)
    for m = 1:M
        divisor -= sims[m].ρ*weights[m]*sims[m].A
        vec_sum += weights[m]*sims[m].δ
    end
    return divisor \ vec_sum
end

function get_submeasure(μ::M, index::AbstractVector{<:Integer}) where M<:AbstractInvariantMeasure
    if index == [0]
        new_μ = μ
    else    
        new_symmetries = μ.symmetries
        # new_suppmeasure = μ.suppmeasure
        for m = index[end:-1:1]
            # new_suppmeasure *= μ.supp.ifs[m].ρ
            new_symmetries = simcompsymmetries(μ.supp.ifs[m], new_symmetries)
        end
        # new_suppmeasure = new_suppmeasure^μ.supp.d
        new_supp = get_subattractor(μ.supp, index)
        new_suppmeasure = μ.suppmeasure * (new_supp.diam / μ.supp.diam) ^ μ.supp.d
        new_barycentre = get_barycentre(new_supp.ifs, μ.weights)
        new_μ = M(new_supp, new_barycentre, new_suppmeasure, μ.weights, new_symmetries)
    end
    return new_μ
end

# is there much argument for defining the submeasure?
# The idea is that it points to the 'parent' anyway,
# and this is memory-cheap because the parent is an immutable struct.
# But V is an immutable struct anyway and μ is necessarilty different in each submeasure.
# So the only issue lies with A.
# dimH and n are the only major duplicates, and these are both tiny in size.
# I also assume that if I write tte constructor carefully, there will be no extra allocations.

# overload the indexing function, so we can get neat vector subscripts
Base.getindex(μ::AbstractInvariantMeasure, inds...) = get_submeasure(μ, [i for i in inds])
Base.getindex(μ::AbstractInvariantMeasure, inds::AbstractVector{<:Integer}) = get_submeasure(μ, inds)
Base.getindex(Γ::AbstractAttractor, inds...) = get_subattractor(Γ, [i for i in inds])
Base.getindex(Γ::AbstractAttractor, inds::AbstractVector{<:Integer}) = get_subattractor(Γ, inds)