function get_bary_weights(μ::AbstractInvariantMeasure{A}, ℓmax::Integer
                        ) where {A<:HomogenousAttractor}
    w = μ.suppmeasure*copy(μ.weights)
    for _ in 2:ℓmax
        w = kron(μ.weights, w)
    end
    return w
end

get_bary_weights(μ::HausdorffMeasure, ℓmax::Integer) = 
    fill(μ.suppmeasure * μ.supp.ρ^(ℓmax*μ.supp.d), length(μ.supp.ifs)^ℓmax)

function barycentre_quadrule( μ::AbstractInvariantMeasure{A}, h::Real
                                ) where {T, R, A<:HomogenousAttractor{T, R}}

    @assert h>0 "Quadrature parameter (second input) must be positive."
    ℓmax = max(ceil(Int64, log(h / μ.supp.diam) / log(μ.supp.ρ)), 0)
    M = length(μ.supp.ifs)
    N = M^ℓmax
    # the above line is the only one which needs modifying for more general measures
    x = Vector{T}(undef, N)
    x[1] = μ.barycentre
    @inbounds for ℓ ∈ 1:ℓmax
        @views x[1:(M^ℓ)] .= μ.supp(x[1:(M^(ℓ-1))])
    end
    return x, get_bary_weights(μ, ℓmax)
end

function barycentre_quadrule(μ₁, μ₂, h)
    x1, w1 = barycentre_quadrule(μ₁, h)
    x2, w2 = barycentre_quadrule(μ₂, h)
    
    return combine_quadrules(x1, w1, x2, w2)
end

barycentre_quadrule(Γ₁::AbstractAttractor, h::Real) =
    barycentre_quadrule(HausdorffMeasure(Γ₁), h::Real)

barycentre_quadrule(Γ₁::AbstractAttractor, Γ₂::AbstractAttractor, h::Real) =
    barycentre_quadrule(HausdorffMeasure(Γ₁), HausdorffMeasure(Γ₂), h::Real)