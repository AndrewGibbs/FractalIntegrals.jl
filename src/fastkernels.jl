function vecdist(x::AbstractVector{<:AbstractVector{T}},
                    y::AbstractVector{<:AbstractVector{T}}
                    ) where T<:Number
    r = zeros(T,size(x))
    n = length(x[1])
    @inbounds @simd for i in eachindex(x)
        @inbounds for j=1:n
            r[i] += (x[i][j]-y[i][j])^2
        end
    end
    return sqrt.(r)
end

function vecdist( x::AbstractVector{T},
                    y::AbstractVector{T}
                    ) where T<:Number
    r = zeros(T,size(x))
    @inbounds @simd for i in eachindex(x)
            r[i] += x[i]-y[i]
        end
    return abs.(r)
end

# likely won't use the following function, but it will be emulated inside quadrature routine
radkernel(  x::AbstractVector,
                y::AbstractVector,
                w::AbstractVector,
                fᵣ::Function
            ) =  @inbounds dot(fᵣ.(vecdist(x,y)),w')

# Laplace/energy kernels
energykernel(s::Number, r::Number) = s==0 ? log(r) : r^-s
energykernel(s::Number, x, y) = energykernel.(s, vecdist(x,y))

# Helmholtz kernels
helmholtzkernel2d(k::Number, x, y) = (im/4)*hankelh1.(0,k.*vecdist(x,y))
function helmholtzkernel3d(k::Number, x, y)
    r = vecdist(x,y)
    return exp.(im*k*r)./(4π*r)
end

# Lipschitz parts of Helmholtz kernels
function helmholtzkernel2d_lipschitzpart(k::Number, x, y)
    x_approx_y = isapprox.(x, y, atol=100*eps(eltype(x)))
    Φ = zeros(ComplexF64, length(x))
    @views Φ[.~x_approx_y] .= helmholtzkernel2d(k, x[.~x_approx_y], y[.~x_approx_y]) .+
                                1/(2π)*energykernel(0.0, x[.~x_approx_y], y[.~x_approx_y])
    Φ[x_approx_y] .= im/4 -1/(2π)*(0.577215664901532 + log(k/2))
    return Φ
end

function helmholtzkernel3d_lipschitzpart(k::Number, x, y)
    x_approx_y = isapprox.(x, y, atol=100*eps(eltype(x)))
    Φ = zeros(ComplexF64, length(x))
    @views r = vecdist(x[.~x_approx_y], y[.~x_approx_y])
    Φ[.~x_approx_y] .= expm1.(im*k*r)./(4π*r)
    Φ[x_approx_y] .= im*k/(4π)
    return Φ
end