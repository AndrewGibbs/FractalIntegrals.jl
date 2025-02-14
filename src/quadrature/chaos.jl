
"""
    x, w = chaos_quadrule(μ::AbstractInvariantMeasure, n::Int)
    x, w = chaos_quadrule(Γ::AbstractAttractor, n::Int)
Randomly allocates quadrature nodes
based on the probability weights of the invariant measure `μ`.
Quadrature weights are uniformly `1/n`

If an attractor Γ is provided, Hausdorff measure is assumed.

This is based on the method introduced in:
Forte, B., Mendivil, F., Vrscay, E.:
"Chaos games for iterated function systems with grey level
maps", 1998.
"""
function chaos_quadrule(
    μ::AbstractInvariantMeasure{<:Any,M},
    numpts::Integer;
    x₀ = get_barycentre(μ),
) where {M}

    # initialise output nodes
    x = Vector{eltype(μ)}(undef, numpts)

    # initial guess, seeds the random algorithm
    x_prev = x₀

    # cumulative sum of probability weights is used to draw random samples
    μ_weights_cum = cumsum(μ.weights)
    for n = 1:numpts
        # choose a random Similarity
        τ = minimum((1:M)[rand().<μ_weights_cum])
        x[n] = μ.supp.ifs[τ](x_prev)

        # log value to map next iteration
        x_prev = x[n]
    end

    # return nodes and equivalued weights
    return x, fill(1 / numpts, numpts)
end

# default to Hausdorff dimension if only an attractor is given
@hausdorffdefault chaos_quadrule

# define for higher dimensional integrals
@tensorquad chaos_quadrule
