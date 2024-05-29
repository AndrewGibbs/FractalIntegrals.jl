# ---------------------- Measures ---------------------------------------------------------#
# abstract type AbstractInvariantMeasure{T<:Real, B, V<:AbstractVector{T}, A<:AbstractAttractor} end

abstract type AbstractInvariantMeasure{A <:AbstractAttractor} end

# AbstractInvariantMeasure{<:AbstractAttractor{T, R}} where T<: Real
# questions: 
    # Is it worth specifying T, R above?
    # Should I define an AbstractSymmetryGroup struct? This would be analagous to AbstractAttractor

struct HausdorffMeasure{T,
                        R <: Real,
                        A <: AbstractAttractor{T, R},
                        V <: AbstractVector{R},
                        G <: AbstractArray{<:AbstractInvariantMap{T}}
                        } <: AbstractInvariantMeasure{A}
    supp::A
    barycentre::T
    suppmeasure::R
    weights::V
    symmetries::G
end

# note the situation with symmetries - they need not be the same for an attractor and a measure.
# But for HausdorffMeasure, the symmetries should be automatically inhereted from the attractor.

struct InvariantMeasure{T,
                        R <: Real,
                        A <: AbstractAttractor{T, R},
                        V <: AbstractVector{T},
                        G <: AbstractArray{<:AbstractInvariantMap{T}}
                        } <: AbstractInvariantMeasure{A}
    supp::A
    barycentre::T
    suppmeasure::R
    weights::V
    symmetries::G
end


InvariantMeasure(Γ::AbstractAttractor,
                weights::AbstractVector{<:Real};
                suppmeasure = 1.0,
                symmetries = trivialgroup(Γ.n)
                ) = 
    InvariantMeasure(Γ,
                    get_barycentre(Γ.ifs, weights),
                    suppmeasure,
                    weights,
                    symmetries)

InvariantMeasure(Γ::AbstractAttractor; vargs...) = HausdorffMeasure(Γ; vargs...)

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

get_hausdorff_weights(Γ::AbstractAttractor) = [s.ρ^Γ.d for s in Γ.ifs]
    
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
        new_suppmeasure = μ.suppmeasure * prod(μ.weights[index])#(new_supp.diam / μ.supp.diam) ^ μ.supp.d
        new_barycentre = get_barycentre(new_supp.ifs, μ.weights)
        new_μ = M(new_supp, new_barycentre, new_suppmeasure, μ.weights, new_symmetries)
    end
    return new_μ
end

# make natural conversions of attractors to Hausdorff measures
HausdorffMeasure(Γ::AbstractAttractor{T, R}) where {T, R} =
    HausdorffMeasure(   Γ,
                        get_barycentre(Γ.ifs, get_hausdorff_weights(Γ)),
                        one(R),
                        get_hausdorff_weights(Γ),
                        Γ.symmetries)

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