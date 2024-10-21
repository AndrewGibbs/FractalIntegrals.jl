# general function for appliying a similarity to an attractor
function (η::AbstractSimilarity)(Γ::A) where A<:AbstractAttractor
    @assert eltype(Γ.ifs) == typeof(η) "Composition similarities must be of same type."
    return A(
            [η ∘ s ∘ η^(-1) for s in Γ.ifs], # new IFS
            η.ρ * Γ.diam, # diameter scaled down
            Γ.d, # dimension unchanged
            Γ.connectedness, # connectedness is unchanged
            [η ∘ s ∘ η^(-1) for s in Γ.symmetries] # symmetries
            )
end

function (η::AbstractSimilarity)(Γ::A) where A<:AbstractHomogenousAttractor
    @assert eltype(Γ.ifs) == typeof(η) "Composition similarities must be of same type."
    return A(
            [η ∘ s ∘ η^(-1) for s in Γ.ifs], # new IFS
            η.ρ * Γ.diam, # diameter scaled down
            Γ.d, # dimension unchanged
            Γ.connectedness, # connectedness is unchanged
            [η ∘ s ∘ η^(-1) for s in Γ.symmetries], # symmetries
            Γ.ρ
            )
end

translating_similarity(x::T) where {T<:Number} =
    OneDimensionalSimilarity{T}(one(T), x, one(T), one(T))

translating_similarity(v::SVector{N, T}) where {N, T} =
    Similarity{N, T}(one(T), v, SMatrix{N, N}(Matrix{T}(IdMat(N))), SMatrix{N, N}(Matrix{T}(IdMat(N))))

# translation by vector
function Base.:+(Γ::AbstractAttractor{N, <:Any, T}, v::AbstractVector{<:Real}) where {N, T}
    @assert N==length(v) "vector must match dimension of attractor"
    η = translating_similarity(SVector{N, T}(v))
    return η(Γ)
end

# translation by scalar
function Base.:+(Γ::AbstractAttractor, x::Real)
    @assert typeof(Γ) <: OneDimensionalAttractorUnion "Attractor must consist of OneDimensional similarities to add scalar"
    η = translating_similarity(x)
    return η(Γ)
end

# commuting case
Base.:+(x, Γ::AbstractAttractor) = Γ + x

function embed_vec(v::SVector{N, T}, n::Integer) where {N, T}
    @assert n ≥ N "must be increasing dimension"
    new_v = zeros(T, n)
    new_v[1:N] .= v
    return SVector{n, T}(new_v)
end

embed_vec(v::T, n::Integer) where T<:Number = embed_vec(SVector{1, T}(v), n)

function embed_mat(m::SMatrix{N, N, T}, n::Integer) where {N, T}
    @assert n ≥ N "must be increasing dimension"
    new_m = Matrix{T}(IdMat(n))
    new_m[1:N,1:N] .= m
    return SMatrix{n, n, T}(new_m)
end

embed_mat(m::T, n::Integer) where T<:Number = embed_mat(SMatrix{1, 1, T}(m), n)

embed(s::AbstractSimilarity, n::Integer) = Similarity(s.ρ,
                                                embed_vec(s.δ, n),
                                                embed_mat(s.A, n),
                                                s.ρ*embed_mat(s.A, n))

# embedding
# function embed(Γ::AbstractAttractor, n::Integer)
#     @assert n ≥ get_ambient_dimension(Γ) "2nd arg must be greater than ambeint fractal dimension"
#     return Attractor(
#             [embed(s, n) for s in Γ.ifs], # new IFS
#             Γ.diam, # diameter scaled down
#             Γ.d, # dimension unchanged
#             Γ.connectedness, # connectedness is unchanged
#             [embed(s, n) for s in Γ.symmetries] # symmetries
#             )
# end

# embedding
function embed(Γ::AbstractAttractor, n::Integer)
    @assert n ≥ get_ambient_dimension(Γ) "2nd arg must be greater than ambeint fractal dimension"
    return Attractor(
            [embed(s, n) for s in Γ.ifs], # new IFS
            diam = Γ.diam, # diameter scaled down
            d = Γ.d, # dimension unchanged
            connectedness = Γ.connectedness, # connectedness is unchanged
            symmetries = [embed(s, n) for s in Γ.symmetries], # symmetries
            )
end

# the types on the left and right aren't the same here!
embed(μ::InvariantMeasure, n::Integer) = InvariantMeasure(embed(μ.supp,n), μ.suppmeasure, μ.weights)
embed(μ::HausdorffMeasure, n::Integer) = HausdorffMeasure(embed(μ.supp,n), μ.suppmeasure, μ.weights)


function embed_into_same_dimension(μple::Tuple{Vararg{AbstractInvariantMeasure}})
    amb_dims = get_ambient_dimension.(μple)
    if fill(amb_dims[1]) == amb_dims
        return μple
    else
        n = maximum(amb_dims)
        return embed.(μple, n)# for μ in μple)
    end
end
# c
# //// /  b n /56trrrffrll;.v“““““““““““““““vvvvvvvvvvv                              fd  hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhnjjjkju