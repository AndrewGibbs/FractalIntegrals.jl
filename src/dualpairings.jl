abstract type AbstractInnerProduct end

struct BarycentreHomogInnerProduct{X<:AbstractVector,
                                    W<:AbstractVector,
                                    S<:AbstractVector,
                                    T<:Real,
                                    I<:AbstractVector{<:T}
                                    } <: AbstractInnerProduct
    x0 :: X
    w :: W
    singular_indices :: S
    singular_integrals :: I
    # log_adjustments :: S
    hQ :: T
end

# outer constructor
function BarycentreHomogInnerProduct(  s::Number,
                                        h_mesh::Real,
                                        Vₕ::FractalBasis{<:HausdorffMeasure{_,_,_,Attr},_},
                                        h_quad::Real
                                        ) where Attr<:HomogenousAttractor

    μ = Vₕ.measure
    Γ = Vₕ.measure.supp

    # first, prepare data for singular integrals
    h_high_scale = Γ.diameter * hQ / h_mesh

    A, B, singular_indices, R, log_adjustments = construct_singularity_matrix(Vₕ.measure, s)
    r = zeros(length(R))
    for n = eachindex(r)
        (m,m_) = R[n]
        x,y,w = barycentre_rule(μ[m], μ[m_], h_high_scale)
        r[n] = w'*energykernel(s,x,y)
    end
    prepared_singular_vals = A\(B*r + log_adjustments) # vector of 'singular values'

    # second, prepare data for smooth integrals
    x, w = barycentre_quadrule(Vₕ[1].measure, h_quad)
    X = x .- Vₕ[1].measure.barycentre

    # create instance, containing everything needed to evaluate dual pairings
    return BarycentreUniformInnerProduct(X, w, singular_indices, prepared_singular_vals, h_quad)
end

function innerproduct(ip::BarycentreHomogInnerProduct, f::Function, ϕ::P0BasisElement)
    x = ip.x0 .+ ϕ.measure.barycentre
    return ϕ.normalisation * dot(f.(x), ip.w')
end

function sesquilinearform(  ip::BarycentreHomogInnerProduct,
                            sio::SeparableIntegralOperator,
                            ϕ::P0BasisElement,
                            ψ::P0BasisElement)
    # first determine if singular or not

    singular_slf = false
    similar_index = 0
    ρ = 0.0

    # first do cheap test - if elements are touching,
    # distance of barycentres must be more than combined diameter:
    if norm(ϕ.barycentre - ψ.barycentre) < (ϕ.diam + ψ.diam)
        # passed initial cheap test for being singular
        # assume that singularity is 'similar' to canonical singular integral
        singular_slf, ρ, similar_index = 
            check_for_similar_integrals(sio.measure,
                                    ip.singular_indices, 
                                    ϕ.vindex,  ψ.vindex,
                                    sio.measure.symmetries, sio.measure.symmetries,
                                    true
                                    )
    end

    # get the tensor product quadrature
    x, y, w = combine_quadrules(ip.x0 .+ ϕ.barycentre, ip.w,
                    ip.x0 .+ ψ.barycentre, ip.w)

    if singular_slf
        scale_adjust, pₙ, pₘ = similar_scaler(ρ,
                                        sio.s,
                                        ip.singular_indices[similar_index][1],
                                        ip.singular_indices[similar_index][2],
                                        ϕ.vindex, ψ.vindex,
                                        ϕ.weights, ψ.weights)

        # compute value of singular integral, as sum of singular + Lipschitz
        I = sio.singularconst * ip.singular_integrals[similar_index] * scale_adjust 
            + dot(sio.lipschitzpart(x,y), conj(w))

        # account for additive log terms, if required
        if sio.s == 0
            I += sio.singularconst * sio.measure.selfmeasure^2 * log(1/ρ) * pₙ * pₘ
        end
    else
        # compute value of smooth integral
        I = dot(sio.kernel(x,y), conj(w))
    end

    # normalise output
    return I * ϕ.normalisation * conj(ψ.normalisation)
end