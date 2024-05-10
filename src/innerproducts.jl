abstract type AbstractInnerProduct end

struct InnerProduct{X<:AbstractVector,
                    W<:AbstractVector,
                    S<:AbstractVector,
                    I<:AbstractVector{<:Number},
                    K<:FractalOperator
                    } <: AbstractInnerProduct
    
    sio :: K
    x :: X
    w :: W
    singular_indices :: S
    singular_integrals :: I
end

# struct HomogInnerProduct{X<:AbstractVector,
#                     W<:Real,
#                     S<:AbstractVector,
#                     # T<:Real,
#                     I<:AbstractVector{<:Number},
#                     K<:FractalOperator
#                     } <: AbstractInnerProduct

#     sio :: K
#     x :: X
#     w :: W
#     singular_indices :: S
#     singular_integrals :: I
# end

function getsingularinfo(μ, s, X, W)
    A, B, singular_indices, R, log_adjustments =
        construct_singularity_matrix(μ, μ, s)
    r = Vector{Float64}(undef,length(R))
    for n = eachindex(r)
        (m,m_) = R[n]
        Xn, Wxn = mapquadrule(μ, m, X, W)
        Yn, Wyn = mapquadrule(μ, m_, X, W)
        x, y, w = combine_quadrules(Xn, Wxn, Yn, Wyn)
        r[n] = w'*energykernel(s, x, y)
    end
    prepared_singular_vals = A\(B*r + log_adjustments) # vector of 'singular values'
    return prepared_singular_vals, singular_indices
end

# outer constructor
# function InnerProduct(  sio::AbstractSingularIntegralOperator,
#                                         Vₕ::FractalBasis{<:HausdorffMeasure{T1,T2,T3,Attr},T4},
#                                         X::AbstractArray,
#                                         W::AbstractArray
#                                         ) where {T1,T2,T3,Attr<:HomogenousAttractor,T4}
function InnerProduct(  sio::AbstractSingularIntegralOperator,
                        Vₕ::FractalBasis,
                        X::AbstractArray,
                        W::AbstractArray)

    μ = Vₕ.measure

    prepared_singular_vals, singular_indices = getsingularinfo(μ, sio.s, X, W)

    # second, prepare data for smooth integrals

    # old homogeneous special case:
    # x, w =  mapquadrule(μ, Vₕ[1].vindex, X, W)
    # X = [xⱼ - Vₕ[1].measure.barycentre for xⱼ in x]
    
    # second, prepare data for smooth integrals
    allquads = Vector{typeof(X)}(undef,length(Vₕ))
    allweights = Vector{typeof(W)}(undef,length(Vₕ))
    for (n, ϕₙ) in enumerate(Vₕ)
        allquads[n], allweights[n] = mapquadrule(μ, ϕₙ.vindex, X, W)
    end

    # create instance, containing everything needed to evaluate dual pairings
    return InnerProduct(sio, allquads, allweights, singular_indices, prepared_singular_vals)
end

# function innerproduct(ip::HomogInnerProduct, f::Function, ψ::P0BasisElement)
#     x = [xⱼ + ψ.measure.barycentre for xⱼ in ip.x0]
#     return conj(ψ.normalisation) * dot(conj(ip.w), f.(x))
# end

innerproduct(ip::InnerProduct, f::Function, ψ::P0BasisElement) = 
    conj(ψ.normalisation) * dot(conj(ip.w[ψ.index]), f.(x[ψ.index]))

function sesquilinearform(  ip::InnerProduct,
                            ϕ::P0BasisElement,
                            ψ::P0BasisElement)
    # first determine if singular or not

    singular_slf = false
    similar_index = 0
    ρ = 0.0

    # first do cheap test - if elements are touching,
    # distance of barycentres must be more than combined diameter:
    if norm(ϕ.measure.barycentre - ψ.measure.barycentre) <
        (ϕ.measure.supp.diam + ψ.measure.supp.diam)
        # passed initial cheap test for being singular
        # assume that singularity is 'similar' to canonical singular integral
        singular_slf, ρ, similar_index = 
            check_for_similar_integrals(ip.sio.measure.supp,
                                    ip.singular_indices, 
                                    ϕ.vindex,
                                    ψ.vindex,
                                    ip.sio.measure.symmetries,
                                    ip.sio.measure.symmetries,
                                    true
                                    )
    end

    # get the tensor product quadrature

    # previous homog special case:
    # x, y, w = combine_quadrules([xⱼ + ϕ.measure.barycentre for xⱼ in ip.x0],
    #                             ip.w[1],
    #                             [xⱼ + ψ.measure.barycentre for xⱼ in ip.x0],
    #                             ip.w[1])

    x, y, w = combine_quadrules(ip.x[ϕ.index],
                                ip.w[ϕ.index],
                                ip.x[ψ.index],
                                ip.w[ψ.index])

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