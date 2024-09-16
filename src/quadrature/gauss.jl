
"""
    x,w = gauss_quadrule(Γ::AbstractInvariantMeasure{<:AbstractAttractor{<:Real, <:Real}},
    N::Integer)

Returns N Gaussian weights w ∈ Rᴺ and nodes x ∈ Rᴺˣᴺ.
Here Γ must be an SelfSimilarFractal in one spatial dimension.
N is the order of the Gauss rule, i.e. number of weights and nodes.
"""
function gauss_quadrule(μ::AbstractInvariantMeasure{AmdDim,<:Any,<:Any,<:AbstractAttractor},
                        N::Integer) where AmdDim
    @assert AmdDim == 1 "Attractor must be compact subset of real line"
    J = getjacobimatrix(μ, N-1)
    vv = real.(eigvecs(J))
    x = real.(eigvals(J))
    w = vv[1,:].^2
    return x,w
end

gauss_quadrule(γ::AbstractAttractor, N::Integer) = gauss_quadrule(HausdorffMeasure(γ),N)