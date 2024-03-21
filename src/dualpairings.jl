abstract type AbstractInnerProduct end

struct BarycentreHomogInnerProduct{X<:AbstractVector,
                                    W<:AbstractVector,
                                    S<:AbstractVector,
                                    T<:Real,
                                    I<:AbstractVector{<:T},
                                    K<:FractalOperator
                                    } <: AbstractInnerProduct
    
    sio :: K
    x0 :: X
    w :: W
    singular_indices :: S
    singular_integrals :: I
    # log_adjustments :: S
    hQ :: T
end

# outer constructor
function BarycentreHomogInnerProduct(  sio::AbstractSingularIntegralOperator,
                                        h_mesh::Real,
                                        Vₕ::FractalBasis{<:HausdorffMeasure{T1,T2,T3,Attr},T4},
                                        hQ::Real
                                        ) where {T1,T2,T3,Attr<:HomogenousAttractor,T4}

    μ = Vₕ.measure
    Γ = Vₕ.measure.supp

    # first, prepare data for singular integrals
    h_high_scale = Γ.diam * hQ / h_mesh

    A, B, singular_indices, R, log_adjustments =
        construct_singularity_matrix(Vₕ.measure, Vₕ.measure, sio.s)
    r = zeros(length(R))
    for n = eachindex(r)
        (m,m_) = R[n]
        x,y,w = barycentre_quadrule(μ[m], μ[m_], h_high_scale)
        r[n] = w'*energykernel(sio.s, x, y)
    end
    prepared_singular_vals = A\(B*r + log_adjustments) # vector of 'singular values'

    # second, prepare data for smooth integrals
    x, w = barycentre_quadrule(Vₕ[1].measure, hQ)
    X = x .- Vₕ[1].measure.barycentre

    # create instance, containing everything needed to evaluate dual pairings
    return BarycentreHomogInnerProduct(sio, X, w, singular_indices, prepared_singular_vals, hQ)
end

function innerproduct(ip::BarycentreHomogInnerProduct, f::Function, ϕ::P0BasisElement)
    x = ip.x0 .+ ϕ.measure.barycentre
    return conj(ϕ.normalisation) * dot(f.(x), conj(ip.w))
end

function sesquilinearform(  ip::BarycentreHomogInnerProduct,
                            ϕ::P0BasisElement,
                            ψ::P0BasisElement)
    # first determine if singular or not

    singular_slf = false
    similar_index = 0
    ρ = 0.0

    # first do cheap test - if elements are touching,
    # distance of barycentres must be more than combined diameter:
    if norm(ϕ.measure.barycentre - ψ.measure.barycentre) < (ϕ.measure.supp.diam + ψ.measure.supp.diam)
        # passed initial cheap test for being singular
        # assume that singularity is 'similar' to canonical singular integral
        singular_slf, ρ, similar_index = 
            check_for_similar_integrals(ip.sio.measure.supp,
                                    ip.singular_indices, 
                                    ϕ.vindex,  ψ.vindex,
                                    ip.sio.measure.symmetries,
                                    ip.sio.measure.symmetries,
                                    true
                                    )
    end

    # get the tensor product quadrature
    x, y, w = combine_quadrules(ip.x0 .+ ϕ.measure.barycentre, ip.w[1],
                                ip.x0 .+ ψ.measure.barycentre, ip.w[1])

    if singular_slf
        scale_adjust, pₙ, pₘ = similar_scaler(ρ,
                                        ip.sio.s,
                                        ip.singular_indices[similar_index][1],
                                        ip.singular_indices[similar_index][2],
                                        ϕ.vindex, ψ.vindex,
                                        ϕ.measure.weights, ψ.measure.weights)

        # compute value of singular integral, as sum of singular + Lipschitz
        I = ip.sio.singularconst * ip.singular_integrals[similar_index] * scale_adjust 
            + dot(ip.sio.lipschitzpart(x,y), conj(w))

        # account for additive log terms, if required
        if ip.sio.s == 0
            I += ip.sio.singularconst * ip.sio.measure.suppmeasure^2 * log(1/ρ) * pₙ * pₘ
        end
    else
        # compute value of smooth integral
        I = dot(ip.sio.kernel(x,y), conj(w))
    end

    # normalise output
    return I * ϕ.normalisation * conj(ψ.normalisation)
end