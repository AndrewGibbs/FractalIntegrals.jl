# ---------------------- Measures ---------------------------------------------------------#
# abstract type AbstractInvariantMeasure{T<:Real, B, V<:AbstractVector{T}, A<:AbstractAttractor} end

abstract type AbstractInvariantMeasure{N, M, T<:Real, A <:AbstractAttractor} end

"""
    InvariantMeasure{N, M, T, A <: AbstractAttractor{N, M, T}}
A general class of measure which can be defined supported on an attractor.
Defined by a 'support', which is an attractor, and a vector of probability weights,
invariant measures satisfy a 'balance' condition.

# Fields
- `supp`: The attractor on which the measure is supported
- `suppmeasure`: The measure of the support
- `weights`: A vector of probability weights
"""
struct InvariantMeasure{N, M, T, A <: AbstractAttractor{N, M, T}
                        } <: AbstractInvariantMeasure{N, M, T, A}
    supp::A
    suppmeasure::T
    weights::SVector{M, T}
end

"""
    HausdorffMeasure{N, M, T, A <: AbstractAttractor{N, M, T}}
A special case of `InvariantMeasure` where
the contraction factors ρ_m and the probability weights p_m satisfy
ρ_m^d=p_m, where d is the Hausdorff dimension of the attractor.
"""
struct HausdorffMeasure{N, M, T, A <: AbstractAttractor{N, M, T}
                        } <: AbstractInvariantMeasure{N, M, T, A}
    supp::A
    suppmeasure::T
    weights::SVector{M, T}
end

# ------------ finding the barycentre (could be moved to barycentre rule?) ------------
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

get_barycentre(μ::AbstractInvariantMeasure) = get_barycentre(μ.supp.ifs, μ.weights)

get_hausdorff_weights(Γ::AbstractAttractor) = [s.ρ^Γ.d for s in Γ.ifs]

# ---------------- outer constructors ----------------------------------


InvariantMeasure(Γ::AbstractAttractor{N, M, T},
                weights::AbstractVector{<:Real};
                suppmeasure = 1.0,
                ) where {N, M, T} = 
    InvariantMeasure(Γ,
                    T(suppmeasure),
                    SVector{M}(T.(weights))
                    )

InvariantMeasure(Γ::AbstractAttractor; vargs...) = HausdorffMeasure(Γ; vargs...)

# make natural conversions of attractors to Hausdorff measures
HausdorffMeasure(Γ::AbstractAttractor{N, M, T}; suppmeasure = one(T)
                ) where {N, M, T} =
    HausdorffMeasure(   Γ,
                        T(suppmeasure),#get_barycentre(Γ.ifs, get_hausdorff_weights(Γ)),
                        SVector{M}(T.(get_hausdorff_weights(Γ)))
                    )

# ------------------------ sub-measure ------------------------------------- #
    
function get_submeasure(μ::M, index::AbstractVector{<:Integer}) where M<:AbstractInvariantMeasure
    if index == [0]
        new_μ = μ
    else    
        new_supp = get_subattractor(μ.supp, index)
        new_suppmeasure = μ.suppmeasure * prod(μ.weights[index])
        new_μ = M(new_supp, new_suppmeasure, μ.weights)
    end
    return new_μ
end

# overload the indexing function, so we can get neat vector subscripts
Base.getindex(μ::AbstractInvariantMeasure, inds::Vararg{<:Integer}) = get_submeasure(μ, [i for i in inds])
Base.getindex(μ::AbstractInvariantMeasure, inds::Vector{<:Integer}) = get_submeasure(μ, inds)
Base.getindex(Γ::AbstractAttractor, inds::Vararg{<:Integer}) = get_subattractor(Γ, [i for i in inds])
Base.getindex(Γ::AbstractAttractor, inds::Vector{<:Integer}) = get_subattractor(Γ, inds)

# there will be lots of cases where we want to default to HausdorffMeasure
macro hausdorffdefault(f)
    quote
        $(esc(f))(Γ::AbstractAttractor, args...; kwargs...) =
            $(esc(f))(HausdorffMeasure(Γ), args...; kwargs...)
    end
end

# It will be useful to extract the type of elements describing the geometry
Base.eltype(μ::AbstractInvariantMeasure) = eltype(μ.supp)

# Some symmetries will be inhereted from attractor, when measure is Hausdorff:
get_symmetries(::InvariantMeasure{N, <:Any, T, <:Any}) where {N, T} = trivialgroup(T, N)
get_symmetries(μ::HausdorffMeasure) = μ.supp.symmetries

get_ambient_dimension(::Union{AbstractAttractor{N}, AbstractInvariantMeasure{N}}) where N = N