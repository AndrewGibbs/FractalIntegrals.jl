# guy here
# https://discourse.julialang.org/t/staticarrays-parametric-types-type-stability/59095
# has same issue as me

# abstract types

abstract type AffineMap{N, T<:Number} end
abstract type AbstractSimilarity{N, T} <: AffineMap{N, T} end

# concrete types

struct Similarity{N, T} <: AbstractSimilarity{N, T}
    ρ :: T
    δ :: SVector{N,T}
    A :: SMatrix{N,N,T}
    ρA :: SMatrix{N,N,T}
end

struct TranslatingSimilarity{N, T} <: AbstractSimilarity{N, T}
    ρ :: T
    δ :: SVector{N,T}
end

struct OneDimensionalSimilarity{N, T} <: AbstractSimilarity{N, T}
    ρ :: T
    δ :: T
    A :: T
    ρA :: T
end

# outer constructors
function Similarity(ρ::T1, δ::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    T = promote_type(T1,T2)
    return TranslatingSimilarity(ρ, SVector{length(δ),T}(δ))
end

function Similarity(ρ::T1, δ::AbstractVector{T2}, A::AbstractMatrix{T3}) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(det(A)) ≈ 1 "Matrix (third argument) must have determinant one"
    N = length(δ)
    @assert (N,N) == size(A) "Matrix dimension must match translation vector"
    T = promote_type(T1,T2,T3)
    return Similarity(ρ, SVector{N,T}(δ), SMatrix{N,N,T}(A), SMatrix{N,N,T}(ρ*A))
end

function Similarity(ρ::T1, δ::T2, r::T3=one(T1)) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(r) ≈ 1 "Reflection (third argument) must have length one"
    T = promote_type(T1,T2,T3)
    return OneDimensionalSimilarity{1,T}(ρ, δ, r, convert(T,ρ*r))
end

# promotion of non-rotating map to rotating map (with identity rotation matrix)

function Base.convert(::Type{Similarity{N, T}}, s::TranslatingSimilarity{N, T}) where {N, T}
    A = Matrix{T}(IdMat(N))
    return Similarity{N, T}(s.ρ, s.δ, A, s.ρ*A)
end

Base.promote_rule(::Type{Similarity{N, T}}, ::Type{TranslatingSimilarity{N, T}}) where {N, T} = Similarity{N, T}

# Similarities as maps, optimised as far as possible

# function (s::Similarity)(x)
#     ρAx = similar(s.δ)
#     mul!(ρAx, s.ρA, x) 
#     return ρAx + s.δ
# end
(s::Similarity)(x)  = s.ρA*x + s.δ
(s::TranslatingSimilarity{N,T})(x) where {N, T} = s.ρ.*SVector{N}(x) .+ s.δ
(s::OneDimensionalSimilarity)(x::Number)  = s.ρA*x + s.δ
sim_map(s::Similarity, x) = s.ρA*x + s.δ

# determine how Similarity appears in the REPL

function Base.show(io::IO, s::TranslatingSimilarity{N,T}) where {N,T}
    str = "x ↦ "*string(round.(s.ρ,digits=2))*"x + "*string(round.(s.δ,digits=2))
    println(typeof(s),':')
    println(io, str)
end

function Base.show(io::IO, s::Similarity)
    str = "x ↦ "*string(round.(s.ρ,digits=2))*string(round.(s.A,digits=2))*"x + "*string(round.(s.δ,digits=2))
    println(typeof(s),':')
    println(io, str)
end