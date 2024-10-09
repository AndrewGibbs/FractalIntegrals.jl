abstract type AbstractInnerProduct end

struct InnerProduct{Q<:AbstractVector{<:QuadStruct},
                    S<:AbstractVector,
                    I<:AbstractVector{<:Number},
                    K<:FractalOperator
                    } <: AbstractInnerProduct
    
    sio :: K
    quadrules :: Q
    singular_indices :: S
    singular_integrals :: I
end

function getsingularinfo(μ, s, q::QuadStruct)
    A, B, singular_indices, R, log_adjustments =
        construct_singularity_matrix(μ, μ, s)
    r = Vector{Float64}(undef, length(R))
    for n = eachindex(r)
        (m,m_) = R[n]
        Xn, Wxn = mapquadrule(μ, m, q.nodes, q.weights)
        Yn, Wyn = mapquadrule(μ, m_, q.nodes, q.weights)
        x, y, w = combine_quadrules(Xn, Wxn, Yn, Wyn)
        r[n] = w'*energykernel(s, x, y)
    end
    prepared_singular_vals = A\(B*r + log_adjustments) # vector of 'singular values'
    return prepared_singular_vals, singular_indices
end

# outer constructor 1
InnerProduct(  sio::AbstractSingularIntegralOperator,
                Vₕ::FractalBasis,
                X::AbstractArray,
                W::AbstractArray) = InnerProduct(sio, Vₕ, QuadStruct(X, W))

# outer constructor 2
function InnerProduct(  sio::AbstractSingularIntegralOperator,
                        Vₕ::FractalBasis,
                        q::QuadStruct)

    μ = Vₕ.measure

    prepared_singular_vals, singular_indices = getsingularinfo(μ, sio.s, q)
    
    # second, prepare data for smooth integrals
    allquads = mapquadrule_to_elements(Vₕ, q)

    # create instance, containing everything needed to evaluate dual pairings
    return InnerProduct(sio, allquads, singular_indices, prepared_singular_vals)
end

innerproduct(ip::InnerProduct, f::Function, ψ::P0BasisElement) = 
    conj(ψ.normalisation) * dot(conj(ip.quadrules[ψ.index].weights), f.(ip.quadrules[ψ.index].nodes))

function sesquilinearform(  ip::InnerProduct,
                            ϕ::P0BasisElement,
                            ψ::P0BasisElement)
    # first determine if singular or not

    singular_slf = false
    similar_index = 0
    ρ = 0.0

    # first do cheap test - if elements are touching,
    # THIS IS NOW QUITE INEFFICIENT AS BOUNDING BALLS ARE COMPUTED O(N^2) TIMES
    if dist⁺(ϕ.measure.supp, ψ.measure.supp) <= 0 #norm(ϕ.measure.barycentre - ψ.measure.barycentre) <
        #(ϕ.measure.supp.diam + ψ.measure.supp.diam)
        # passed initial cheap test for being singular
        # assume that singularity is 'similar' to canonical singular integral
        
        parent_symmetries = get_symmetries(ip.sio.measure)

        singular_slf, ρ, similar_index = 
            check_for_similar_integrals(ip.sio.measure.supp,
                                    ip.singular_indices, 
                                    Vector(ϕ.vindex),
                                    Vector(ψ.vindex),
                                    parent_symmetries,
                                    parent_symmetries,
                                    true
                                    )
    end

    # get the tensor product quadrature

    x, y, w = combine_quadrules(ip.quadrules[ϕ.index].nodes,
                                ip.quadrules[ϕ.index].weights,
                                ip.quadrules[ψ.index].nodes,
                                ip.quadrules[ψ.index].weights)

    if singular_slf
        scale_adjust = similar_scaler(ρ,
                                    ip.sio.s,
                                    ip.singular_indices[similar_index][1],
                                    ip.singular_indices[similar_index][2],
                                    ϕ.vindex, ψ.vindex,
                                    ϕ.measure.weights, ψ.measure.weights)

        # compute value of singular integral, as sum of singular + Lipschitz
        I = ip.sio.singularconst * ip.singular_integrals[similar_index] * scale_adjust +
            dot(conj(w), ip.sio.lipschitzpart(x,y))

        # account for additive log terms, if required
        if ip.sio.s == 0
            pϕ = compose_weights(ϕ.measure.weights, ϕ.vindex)
            pψ = compose_weights(ψ.measure.weights, ψ.vindex)
            I += ip.sio.singularconst * ip.sio.measure.suppmeasure^2 * log(1/ρ) * pϕ * pψ
        end
    else
        # compute value of smooth integral
        I = dot(conj(w), ip.sio.kernel(x,y))
    end

    # normalise output
    return I * ϕ.normalisation * conj(ψ.normalisation)
end