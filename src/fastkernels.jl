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

energykernel(s::Number, r::Number) = s==0 ? log(r) : r^-s
energykernel(s::Number, x, y) = energykernel.(s, vecdist(x,y))