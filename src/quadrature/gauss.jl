
"""
    x, w = gauss_quadrule(μ::AbstractInvariantMeasure, n::Integer)
    x, w = gauss_quadrule(Γ::AbstractAttractor, n::Integer)
Returns `n` Gaussian weights `w ∈ ℝⁿ` and nodes `x ∈ ℝⁿ`.
Here n is the order of the Gauss rule, i.e. number of weights and nodes.

This is based on the algorithm introduced in:
"A stable Stieltjes technique for computing orthogonal polynomials
    and Jacobi matrices associated with a class of singular measures",
Mantica, 1996.
As in this paper, the algorithm is only defined for measures and
attractors in one ambient dimension.
"""
function gauss_quadrule(μ::AbstractInvariantMeasure{AmdDim,<:Any,<:Any,<:AbstractAttractor},
                        numpts::Integer) where AmdDim
    @assert AmdDim == 1 "Attractor must be compact subset of real line"
    J = getjacobimatrix(μ, numpts-1)
    vv = real.(eigvecs(J))
    x = real.(eigvals(J))
    w = vv[1,:].^2
    return x, w
end

@hausdorffdefault gauss_quadrule
# gauss_quadrule(γ::AbstractAttractor, N::Integer) = gauss_quadrule(HausdorffMeasure(γ),N)

@tensorquad gauss_quadrule