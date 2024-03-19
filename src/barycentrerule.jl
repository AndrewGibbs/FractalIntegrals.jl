

function barycentre_quadrule( μ::F, h::H
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
    w = fill(μ.suppmeasure * μ.supp.ρ^(ℓmax*μ.supp.d), N)
    # the above line is the only one which needs modifying for more general measures
    x = zeros(B, N)
    x[1] = μ.barycentre
    @inbounds for ℓ ∈ 1:ℓmax
        @views x[1:(M^ℓ)] .= μ.supp(x[1:(M^(ℓ-1))])
    end
    return x, w
end

# the below function should be in a more generic location.
# also, there should be a fast version for scalar w, i.e. when all weights are equal
function combine_quadrules(x1, w1, x2, w2)
return repeat(x1, inner=length(w2)),
        repeat(x2, outer=length(w1)),
        repeat(w1, inner=length(w2)).*repeat(w2,outer=length(w1))

function combine_quadrules(x1::AbstractArray{<:SVector},
                        w1<:AbstractArray{<:Number},
                        x2::AbstractArray{<:SVector},
                        w2<:AbstractArray{<:Number})
    N1 = length(w1)
    N2 = length(w2)
    N = N1*N2
    X1 = zeros(eltype(x1), N)
    X2 = zeros(eltype(x2), N)
    W1 = zeros(eltype(w1), N)
    W2 = zeros(eltype(w2), N)

    @simd for n1 in 1:N1
        X1[((n1-1)*N2+1):(n1*N2)] .= x1
        W1[((n1-1)*N2+1):(n1*N2)] .= w1
    end

    @simd for n2 in 1:N2
        X2[((n2-1)*N1+1):(n2*N1)] .= x2[n2]
        W2[((n2-1)*N1+1):(n2*N1)] .= w2[n2]
    end

    return X1, X2, W1.*W2
end

function combine_quadrules(x1::AbstractArray{<:SVector},
    w1::Number,
    x2::AbstractArray{<:SVector},
    w2::Number)

    N1 = length(x1)
    N2 = length(x2)
    N = N1*N2
    X1 = zeros(eltype(x1), N)
    X2 = zeros(eltype(x2), N)
    # W1 = fill(eltype(w1), N)
    # W2 = zeros(eltype(w2), N)

    @simd for n1 in 1:N1
        X1[((n1-1)*N2+1):(n1*N2)] .= x1
    # W1[((n1-1)*N2+1):(n1*N2)] .= w1
    end

    @simd for n2 in 1:N2
        X2[((n2-1)*N1+1):(n2*N1)] .= x2[n2]
    # W2[((n2-1)*N1+1):(n2*N1)] .= w2[n2]
    end

    return X1, X2, fill(w1*w2, N)
end

function barycentre_quadrule(μ₁, μ₂, h)
    x1, w1 = barycentre_quadrule(μ₁, h)
    x2, w2 = barycentre_quadrule(μ₂, h)
    
    return combine_quadrules(x1, w1, x2, w2)
end