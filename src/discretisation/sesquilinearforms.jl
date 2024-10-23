(φ::PreQuadP0BasisElement)(f::Function) = φ.quadrule.weights' * (f.(φ.quadrule.nodes) .* conj.(φ.(φ.quadrule.nodes)))

dot(φ::PreQuadP0BasisElement, f::Function) = φ(f)

function dot(ϕ::AbstractP0BasisElement, ψ::AbstractP0BasisElement)
    # start with lengths of vector indices
    ϕ_leng = length(ϕ.vindex)
    ψ_leng = length(ψ.vindex)

    # if fractal depth is the same
    if ϕ_leng == ψ_leng
        # check if vector indices match
        if ϕ.vindex == ψ.vindex
            overlap_measure = ϕ.measure.suppmeasure
        else
            overlap_measure = zero(ϕ.measure.suppmeasure)
        end
    else # if depths are different, one fractal may be subset of the other
        least_deep_basfn = [ϕ, ψ][argmin(ϕ_leng, ψ_leng)]
        most_deep_basfn = [ϕ, ψ][argmax(ϕ_leng, ψ_leng)]
        if most_deep_basfn.vindex[1:length(least_deep_basfn)] == least_deep_basfn
            overlap_measure = least_deep_basfn.measure.suppmeasure
        else
            overlap_measure = zero(ϕ.measure.suppmeasure)
        end
    end

    return overlap_measure * ϕ.normalisation * conj(ψ.normalisation)
end

(φ::PreQuadP0BasisElement)(ψ::PreQuadP0BasisElement) = φ⋅ψ

# machinery for Galerkin matrix construction 

function sesquilinearform(  op::SmoothIntegralOperator,
                            ϕ::PreQuadP0BasisElement,
                            ψ::PreQuadP0BasisElement)
    x, y, w = combine_quadrules(ϕ.quadrule.nodes,
                                ϕ.quadrule.weights,
                                ψ.quadrule.nodes,
                                ψ.quadrule.weights)
    return dot(conj(w), op.kernel(x,y))
end

function sesquilinearform(  sio::AbstractSeparableIntegralOperator,
                            ϕ::PreQuadP0BasisElement,
                            ψ::PreQuadP0BasisElement,
                            singular_integrals::Vector{<:Number},
                            singular_indices::Vector{<:Tuple{Vector{<:Integer}, Vector{<:Integer}}})

    singular_slf = false
    similar_index = 0
    ρ = 0.0

    if dist⁻(ϕ.measure.supp, ψ.measure.supp) <= 0 #norm(ϕ.measure.barycentre - ψ.measure.barycentre) <
        #(ϕ.measure.supp.diam + ψ.measure.supp.diam)
        # passed initial cheap test for being singular
        # assume that singularity is 'similar' to canonical singular integral

        parent_symmetries = get_symmetries(sio.measure)

        singular_slf, ρ, similar_index = 
        check_for_similar_integrals(sio.measure.supp,
                        singular_indices, 
                        Vector(ϕ.vindex),
                        Vector(ψ.vindex),
                        parent_symmetries,
                        parent_symmetries,
                        true
                        )
    end

        # get the tensor product quadrature

        x, y, w = combine_quadrules(ϕ.quadrule.nodes,
                    ϕ.quadrule.weights,
                    ψ.quadrule.nodes,
                    ψ.quadrule.weights)

    if singular_slf
    scale_adjust = similar_scaler(ρ,
                    sio.s,
                    singular_indices[similar_index][1],
                    singular_indices[similar_index][2],
                    ϕ.vindex, ψ.vindex,
                    ϕ.measure.weights, ψ.measure.weights)

    # compute value of singular integral, as sum of singular + Lipschitz
    I = sio.singularconst * singular_integrals[similar_index] * scale_adjust +
    dot(conj(w), sio.lipschitzpart(x,y))

        # account for additive log terms, if required
        if sio.s == 0
            pϕ = compose_weights(ϕ.measure.weights, ϕ.vindex)
            pψ = compose_weights(ψ.measure.weights, ψ.vindex)
            I += sio.singularconst * sio.measure.suppmeasure^2 * log(1/ρ) * pϕ * pψ
        end
    else
    # compute value of smooth integral
        I = dot(conj(w), sio.kernel(x,y))
    end

    # normalise output
    return I * ϕ.normalisation * conj(ψ.normalisation)
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