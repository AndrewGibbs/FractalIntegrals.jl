

function barycentre_quadrature( μ::F, h::H
                                ) where {
                                T,
                                B,
                                V,
                                A<:HomogenousAttractor,
                                H<:Real, 
                                F<:HausdorffMeasure{T,B,V,A}
                                }

    ℓmax = ceil(Int64, log(h / μ.supp.diam) / log(μ.supp.ρ))
    M = length(μ.supp.ifs)
    N = M^ℓmax
    # consider splitting functions here - so container sizes can be determined by compiler
    w = fill(μ.suppmeasure * μ.supp.ρ^(ℓmax*μ.supp.d), N)
    x = zeros(B, N)
    x[1] = μ.barycentre
    @inbounds for ℓ ∈ 1:ℓmax
        @views x[1:(M^ℓ)] .= μ.supp(x[1:(M^(ℓ-1))])
        # ifs_map!(x[1:(M^ℓ)], μ.supp.ifs, copy(x[1:(M^(ℓ-1))]))
    end
    return x, w
end