struct ColPt{V<:VectorIndex,T,S<:Real}
    vindex::V
    node::T
    weight::S
end

struct DiscreteCollocationOperator{
    K<:FractalOperator,
    B<:FractalBasis,
    M<:AbstractMatrix,
    T<:AbstractVector{<:ColPt},
} <: DiscreteFractalOperator
    op::K
    basis::B
    collocation_matrix::M
    collocation_points::T
end

function get_collocation_points(μ::AbstractInvariantMeasure, h_mesh::Real)
    Lₕ = subdivide_indices(μ.supp, h_mesh)
    return [ColPt(𝐦, get_barycentre(μ[𝐦]), μ[𝐦].suppmeasure) for 𝐦 in Lₕ]
end

function get_collocation_points(Vₕ::FractalBasis)
    return [
        ColPt(ϕₕ.vindex, get_barycentre(ϕₕ.measure), ϕₕ.measure.suppmeasure) for ϕₕ in Vₕ
    ]
end

# this is basically a graded quadrature routine. Should be generalised and kept elsewhere.
function compute_col_entry(
    op::IntegralOperator,
    ϕ::FractalBasisElement,
    x::ColPt,
    h_quad::Real;
    min_mesh_width_permitted = 1e-2,
)
    𝓖ₕ, condish_satisfied = grade_mesh_towards_point(
        ϕ.measure.supp,
        x.node;
        h = h_quad,
        min_mesh_width_permitted = min_mesh_width_permitted,
    )

    if dist⁻(ϕ.measure.supp, x.node) <= 0 # node is close to support of basis element
        graded_mesh = [ϕ.measure[𝐤] for 𝐤 in 𝓖ₕ]
        Φₓ(y) = op.kernel(x.node, y)

        # get nodes (barycentres) and weights (measures) of elements.
        inner_quad_nodes = [get_barycentre(γ) for γ in graded_mesh][condish_satisfied]
        inner_quad_weights = [γ.suppmeasure for γ in graded_mesh][condish_satisfied]

        # now apply quadrature rule to kernel and return value
        # return conj(inner_quad_weights)' ⋅ Φₓ.(inner_quad_nodes)
        return conj(inner_quad_weights)' ⋅
               (op.kernel(x.node, inner_quad_nodes) .* ϕ(inner_quad_nodes))
    else
        return conj(ϕ.quadrule.weights)' ⋅
               (op.kernel(x.node, ϕ.quadrule.nodes) .* ϕ(ϕ.quadrule.nodes))
    end
end

# value of ϕ at the collocation point: non-zero when the point's element is inside ϕ's
function compute_col_entry(
    ::IdentityOperator,
    ϕ::FractalBasisElement,
    x::ColPt,
    ::Real;
    varargs...,
)
    𝐦, 𝐧 = Vector(ϕ.vindex), Vector(x.vindex)
    nested = 𝐦 == [0] || (length(𝐦) ≤ length(𝐧) && 𝐧[1:length(𝐦)] == 𝐦)
    return nested ? ϕ.normalisation : zero(ϕ.normalisation)
end

function compute_col_entry(
    S::ScaledOperator,
    ϕ::FractalBasisElement,
    x::ColPt,
    h_quad::Real;
    varargs...,
)
    return S.λ * compute_col_entry(S.operator, ϕ, x, h_quad; varargs...)
end

function compute_col_entry(
    S::SumOperator,
    ϕ::FractalBasisElement,
    x::ColPt,
    h_quad::Real;
    varargs...,
)
    return compute_col_entry(S.operator1, ϕ, x, h_quad; varargs...) +
           compute_col_entry(S.operator2, ϕ, x, h_quad; varargs...)
end

function get_collocation_matrix(
    op::FractalOperator,
    basis::FractalBasis,
    h_col::Real,
    h_quad::Real;
    varargs...,
)
    num_basis_els = length(basis)
    h_col < Inf ? oversample = true : oversample = false

    if oversample
        col_pts = get_collocation_points(op.measure, h_col)
    else
        col_pts = get_collocation_points(basis)
    end
    num_col_pts = length(col_pts)

    @assert num_basis_els <= num_col_pts "Number of collocation points must be no less than basis size"

    col_matrix = Matrix{ComplexF64}(undef, num_col_pts, num_basis_els)

    @sync for m = 1:num_col_pts
        @spawn for n = 1:num_basis_els
            @inbounds col_matrix[m, n] =
                compute_col_entry(op, basis[n], col_pts[m], h_quad; varargs...)
        end
    end

    return col_matrix, col_pts
end

# generic discretisation function
function discretise_collocation(
    K::FractalOperator,
    Vₕ::FractalBasis,
    h_quad;
    h_col = Inf,
    kwargs...,
)
    return DiscreteCollocationOperator(
        K,
        Vₕ,
        get_collocation_matrix(K, Vₕ, h_col, h_quad; kwargs...)...,
    )
end

# truncated SVD algorithm
function pseudo_backslash(A, b; ϵ = sqrt(eps()))
    U, Σ, V = svd(A)

    pseudo_inv_diag = Diagonal(zeros(length(Σ)))
    for n = 1:length(Σ)
        if abs(Σ[n]) / abs(Σ[1]) > ϵ
            pseudo_inv_diag[n, n] = 1 / Σ[n]
        end
    end

    return V * pseudo_inv_diag * U' * b
end

function Base.:\(op::DiscreteCollocationOperator, f::Function; varargs...)
    # get collocation samples of RHS
    f_samples = [f(x.node) for x in op.collocation_points]

    # determine if oversampling has been used
    if length(op.collocation_points) > length(op.basis)
        # if so, use truncated SVD
        coeffs = pseudo_backslash(op.collocation_matrix, f_samples; varargs...)
    elseif length(op.collocation_points) == length(op.basis)
        # otherwise, standard matrix inversion
        coeffs = op.collocation_matrix \ f_samples
    else
        # check for reasonable matrix shape
        error("Fewer collocation points than basis elements")
    end
    return Projection(op.basis, coeffs)
end
