function get_bary_weights(
    μ::AbstractInvariantMeasure{<:Any,<:Any,<:Any,<:HomogenousAttractor},
    ℓmax::Integer,
)
    w = μ.suppmeasure * copy(μ.weights)
    for _ = 2:ℓmax
        w = kron(μ.weights, w)
    end
    return w
end

function get_bary_weights(μ::HausdorffMeasure, ℓmax::Integer)
    return fill(μ.suppmeasure * μ.supp.ρ^(ℓmax * μ.supp.d), length(μ.supp.ifs)^ℓmax)
end

"""
    x, w = barycentre_quadrule(μ::AbstractInvariantMeasure, h::Real)
    x, w = barycentre_quadrule(Γ::AbstractAttractor, h::Real)
Subdivides the support `Γ` of the measure `μ` into self-similar components,
no bigger than diameter `h`,
and allocates quadrature points at the barycentre of each,
with weights corresponding to the respective measure of each subcomonent.

If an attractor Γ is provided, Hausdorff measure is assumed.

This is based on the method introduced in:
"Numerical quadrature for singular integrals on fractals",
A. Gibbs, D. P. Hewett, A. Moiola, 2022.
"""
function barycentre_quadrule(
    μ::AbstractInvariantMeasure{<:Any,M,<:Any,<:HomogenousAttractor},
    h::Real,
) where {M}
    @assert h > 0 "Quadrature parameter (second input) must be positive."
    ℓmax = max(ceil(Int64, log(h / μ.supp.diam) / log(μ.supp.ρ)), 0)
    # M = length(μ.supp.ifs)
    max_num_pts = M^ℓmax
    # the above line is the only one which needs modifying for more general measures
    x = Vector{eltype(μ)}(undef, max_num_pts)
    x[1] = get_barycentre(μ)
    @inbounds for ℓ ∈ 1:ℓmax
        @views x[1:(M^ℓ)] .= μ.supp(x[1:(M^(ℓ-1))])
    end
    return x, get_bary_weights(μ, ℓmax)
end

function get_bary_leaf(𝐦, adj_tree_dict, μ)
    @assert haskey(adj_tree_dict, [0]) "dict must have zero index [0]"

    # simplify local notation
    S = μ.supp.ifs
    w = μ.weights

    # start by looking for the longest [j:end] index string in dict
    start_index = length(𝐦)
    𝐦_in_dict = Int64[]
    loop_break = false

    # find the largest useful sub-word which is already in the dict 
    for i in eachindex(𝐦)
        𝐦_check = 𝐦[i:end]
        if haskey(adj_tree_dict, 𝐦_check)
            start_index = i - 1
            𝐦_in_dict = 𝐦_check
            loop_break = true
            break
        end
    end

    # account for case where dict is effectively empty and needs first (nontrivial) entry
    if !loop_break
        adj_tree_dict[[𝐦[end]]] = ( # index
            S[𝐦[end]](adj_tree_dict[[0]][1]), # nodes
            w[𝐦[end]] * adj_tree_dict[[0]][2], # weights
        )
        start_index -= 1
    end

    # now work down the adjoint tree, filling in words into the dict
    𝐦_in_dict = 𝐦[(start_index+1):end]
    for i = start_index:-1:1
        𝐦_new = 𝐦[i:end]#[i; 𝐦_in_dict]
        adj_tree_dict[𝐦_new] = ( # index
            S[𝐦[i]](adj_tree_dict[𝐦_in_dict][1]), # nodes
            w[𝐦[i]] * adj_tree_dict[𝐦_in_dict][2], # weights
        )
        𝐦_in_dict = 𝐦_new
    end

    # return the weights and nodes for the requested word
    return adj_tree_dict[𝐦]
end

function barycentre_quadrule(
    μ::AbstractInvariantMeasure{N,M,T,A},
    h::Number,
) where {N,M,T,A}
    NodeType = A <: OneDimensionalAttractorUnion ? T : SVector{N,T}
    Lₕ = subdivide_indices(μ.supp, h)

    # initialise tree
    num_nodes = length(Lₕ)
    adj_tree_dict = Dict{Vector{Int64},Tuple{NodeType,T}}()

    # allocate initial barycentre
    adj_tree_dict[[0]] = (get_barycentre(μ), μ.suppmeasure)

    # allocate output data
    nodes = Vector{NodeType}(undef, num_nodes)
    weights = Vector{T}(undef, num_nodes)
    for (n, 𝐦) in enumerate(Lₕ)
        nodes[n], weights[n] = get_bary_leaf(𝐦, adj_tree_dict, μ)
    end

    return nodes, weights
end

function barycentre_quadrule(μ₁, μ₂, h)
    x1, w1 = barycentre_quadrule(μ₁, h)
    x2, w2 = barycentre_quadrule(μ₂, h)

    return combine_quadrules(x1, w1, x2, w2)
end

@hausdorffdefault barycentre_quadrule
# barycentre_quadrule(Γ₁::AbstractAttractor, h::Real) =
#     barycentre_quadrule(HausdorffMeasure(Γ₁), h::Real)

# barycentre_quadrule(Γ₁::AbstractAttractor, Γ₂::AbstractAttractor, h::Real) =
#     barycentre_quadrule(HausdorffMeasure(Γ₁), HausdorffMeasure(Γ₂), h::Real)

@tensorquad barycentre_quadrule
