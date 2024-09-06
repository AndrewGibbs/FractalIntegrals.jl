# ---------------------- Measures ---------------------------------------------------------#
# abstract type AbstractInvariantMeasure{T<:Real, B, V<:AbstractVector{T}, A<:AbstractAttractor} end

abstract type AbstractInvariantMeasure{N, M, T<:Real, A <:AbstractAttractor} end
 
struct InvariantMeasure{N, M, T, A <: AbstractAttractor{N, M, T}
                        } <: AbstractInvariantMeasure{N, M, T, A}
    supp::A
    suppmeasure::T
    weights::SVector{M, T}
end

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
        # new_symmetries = μ.symmetries
        # # new_suppmeasure = μ.suppmeasure
        # for m = index[end:-1:1]
        #     # new_suppmeasure *= μ.supp.ifs[m].ρ
        #     new_symmetries = simcompifs(μ.supp.ifs[m], new_symmetries)
        #     # new_symmetries = [simcomp(μ.supp.ifs[m], ns) for ns in new_symmetries]#simcompsymmetries(μ.supp.ifs[m], new_symmetries)
        # end
        # new_suppmeasure = new_suppmeasure^μ.supp.d
        new_supp = get_subattractor(μ.supp, index)
        new_suppmeasure = μ.suppmeasure * prod(μ.weights[index])#(new_supp.diam / μ.supp.diam) ^ μ.supp.d
        # new_barycentre = get_barycentre(new_supp.ifs, μ.weights)
        # new_μ = M(new_supp, new_barycentre, new_suppmeasure, μ.weights, new_symmetries)
        new_μ = M(new_supp, new_suppmeasure, μ.weights)
    end
    return new_μ
end

# overload the indexing function, so we can get neat vector subscripts
Base.getindex(μ::AbstractInvariantMeasure, inds...) = get_submeasure(μ, [i for i in inds])
Base.getindex(μ::AbstractInvariantMeasure, inds::AbstractVector{<:Integer}) = get_submeasure(μ, inds)
Base.getindex(Γ::AbstractAttractor, inds...) = get_subattractor(Γ, [i for i in inds])
Base.getindex(Γ::AbstractAttractor, inds::AbstractVector{<:Integer}) = get_subattractor(Γ, inds)

# there will be lots of cases where we want to default to HausdorffMeasure

# Use metaprogramming to do this
# macro default_to_hausdorff(funcs...)
#     quote
#         $(Expr(:block, [:(function $func(Γ::AbstractAttractor, varargs...; kwargs...)
#                               $func(HausdorffMeasure(Γ), varargs...; kwargs...)
#                           end) for func in funcs]...))
#     end
# end