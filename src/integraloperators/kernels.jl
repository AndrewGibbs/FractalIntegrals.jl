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

function vecdist(x::AbstractVector{T},
    y::AbstractVector{<:AbstractVector{T}}
    ) where T<:Number
    r = zeros(T,size(y))
    n = length(x)
    @inbounds @simd for i in eachindex(y)
        @inbounds for j=1:n
            r[i] += (x[j]-y[i][j])^2
        end
    end
    return sqrt.(r)
end

# distance commutes
vecdist(x::AbstractVector{<:AbstractVector{T}}, y::AbstractVector{T}
        ) where T<:Number = vecdist(y, x)

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
helmholtzkernel2d(k::Number, x, y) = (im/4)*hankelh1.(0,k*vecdist(x,y))
function helmholtzkernel3d(k::Number, x, y)
    r = vecdist(x,y)
    return exp.(im*k*r)./(4π*r)
end

# Lipschitz parts of Helmholtz kernels
function helmholtzkernel2d_lipschitzpart(k::Number, x, y)
    T = eltype(x[1])
    x_approx_y = isapprox.(x, y, atol=100*eps(T))
    Φ = zeros(Complex{T}, length(x))
    if false ∈ x_approx_y
        @views Φ[.~x_approx_y] .= helmholtzkernel2d(k, x[.~x_approx_y], y[.~x_approx_y]) .+
                                    1/(2π)*energykernel(0.0, x[.~x_approx_y], y[.~x_approx_y])
    end
    if true ∈ x_approx_y
        Φ[x_approx_y] .= im/4 -1/(2π)*(0.577215664901532 + log(k/2))
    end
    return Φ
end

function helmholtzkernel3d_lipschitzpart(k::Number, x, y)
    T = eltype(x[1])
    x_approx_y = isapprox.(x, y, atol=100*eps(T))
    Φ = Vector{Complex{T}}(undef, length(x))
    if false ∈ x_approx_y
        @views r = vecdist(x[.~x_approx_y], y[.~x_approx_y])
        Φ[.~x_approx_y] .= expm1.(im*k*r)./(4π*r)
    end
    if true ∈ x_approx_y
        Φ[x_approx_y] .= im*k/(4π)
    end
    return Φ
end