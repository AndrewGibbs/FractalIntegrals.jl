# guy here
# https://discourse.julialang.org/t/staticarrays-parametric-types-type-stability/59095
# has same issue as me

# abstract types

abstract type AffineMap end
abstract type AbstractSimilarity <: AffineMap end

# concrete types

struct Similarity{T<:Number, V<:AbstractVector{T}, M<:AbstractMatrix{T}} <: AbstractSimilarity
    ρ :: T
    δ :: V
    A :: M
    ρA :: M
end

struct TranslatingSimilarity{T<:Number, V<:AbstractVector{T}} <: AbstractSimilarity
    ρ :: T
    δ :: V
end

struct OneDimensionalSimilarity{T<:Number} <: AbstractSimilarity
    ρ :: T
    δ :: T
    A :: T
    ρA :: T
end

# outer constructors
function Similarity(ρ::T1, δ::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    T = promote_type(T1,T2)
    return TranslatingSimilarity(T(ρ), SVector{length(δ),T}(δ))
end

function Similarity(ρ::T1, δ::AbstractVector{T2}, A::AbstractMatrix{T3}) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(det(A)) ≈ 1 "Matrix (third argument) must have determinant one"
    N = length(δ)
    @assert (N,N) == size(A) "Matrix dimension must match translation vector"
    T = promote_type(T1,T2,T3)
    #{T,SVector{N,T},SMatrix{N,N,T,N^2}(A)}
    return Similarity(T(ρ), SVector{N,T}(δ), SMatrix{N,N,T}(A), SMatrix{N,N,T}(ρ*A))
end

function Similarity(ρ::T1, δ::T2, r::T3=one(T1)) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(r) ≈ 1 "Reflection (third argument) must have length one"
    T = promote_type(T1,T2,T3)
    return OneDimensionalSimilarity(T(ρ), T(δ), T(r), convert(T,ρ*r))
end

# promotion of non-rotating map to rotating map (with identity rotation matrix)

function Base.convert(::Type{Similarity{T,V,M}}, s::TranslatingSimilarity{T,V}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix}
    A = convert(M,IdMat(length(s.δ)))
    return Similarity(s.ρ, s.δ, A, s.ρ*A)
end

# Base.promote_rule(Type{Similarity{T, V, M}}, Type{TranslatingSimilarity{T, V}}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix} = Similarity{T, V, M}
Base.promote_rule(::Type{Similarity{T, V, M}}, ::Type{TranslatingSimilarity{T, V}}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix} = Similarity{T, V, M}

# Similarities as maps, optimised as far as possible

(s::AbstractSimilarity)(x)  = s.ρA*x + s.δ
(s::TranslatingSimilarity)(x) = s.ρ*convert(typeof(s.δ),x) + s.δ

# function ifs_map!(  Sx::AbstractVector, #output
#                     S::AbstractVector{<:AbstractSimilarity},
#                     x::AbstractVector)
#     N = length(x)
#     for m in eachindex(S)
#         # @simd for j in eachindex(x)
#         @views Sx[N*(m-1)+1 : N*m] .= S[m].ρA.*x .+ S[m].δ
#         # end
#     end
# end

# function ifs_map(S::AbstractVector{<:AbstractSimilarity},
#                 x::AbstractVector{T}) where T
#     Sx = zeros(T,length(S)*length(x))
#     ifs_map!(Sx, S, x)
#     return Sx
# end

# can't appear to gain anything by using in-place matvecs
# function sim_map!(s::AbstractSimilarity, x)
#     mul!(x,s.ρA,x)
#     x .+= s.δ
# end

# determine how Similarity appears in the REPL

info_string(s::TranslatingSimilarity) = "x ↦ "*string(round.(s.ρ,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::Similarity) = "x ↦ "*string(round.(s.ρ,digits=2))*string(round.(s.A,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::OneDimensionalSimilarity) = "x ↦ "*string(round(s.ρA,digits=2))*"x + "*string(round(s.δ,digits=2))

Base.show(io::IO, s::AbstractSimilarity)  = print(io,'\n',typeof(s),':','\n', info_string(s))

# affine maps

abstract type AbstractInvariantMap <:AffineMap end
struct IdentityMap <: AbstractInvariantMap end

(I::IdentityMap)(x) = x