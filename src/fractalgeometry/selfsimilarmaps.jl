abstract type AbstractSimilarity{R<:Real, T} <: AffineMap end

# concrete types

struct Similarity{R, T, M<:AbstractMatrix} <: AbstractSimilarity{R, T}
    ρ :: R
    δ :: T
    A :: M
    ρA :: M
end

struct TranslatingSimilarity{R, T} <: AbstractSimilarity{R, T}
    ρ :: R
    δ :: T
    A :: R
    ρA :: R
end

struct OneDimensionalSimilarity{R, T} <: AbstractSimilarity{R, T}
    ρ :: R
    δ :: T
    A :: T
    ρA :: T
end

# outer constructors
function Similarity(ρ::T1, δ::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    T = promote_type(T1, T2)
    return TranslatingSimilarity(T(ρ), SVector{length(δ),T}(δ), one(T), T(ρ))
end

function Similarity(ρ::T1, δ::AbstractVector{T2}, A::AbstractMatrix{T3}
                    ) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(det(A)) ≈ 1 "Matrix (third argument) must have determinant one"
    N = length(δ)
    @assert (N,N) == size(A) "Matrix dimension must match translation vector"
    T = promote_type(T1,T2,T3)
    return Similarity(T(ρ), SVector{N,T}(δ), SMatrix{N,N,T}(A), SMatrix{N,N,T}(ρ*A))
end

function Similarity(ρ::T1, δ::T2, r::T3=one(T1)) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(r) ≈ 1 "Reflection (third argument) must have length one"
    T = promote_type(T1,T2,T3)
    return OneDimensionalSimilarity(T(ρ), T(δ), T(r), convert(T,ρ*r))
end

# promotion of non-rotating map to rotating map (with identity rotation matrix)

function Base.convert(::Type{Similarity{R, T, M}}, s::TranslatingSimilarity{R, T}
                        ) where {R<:Real, T<:AbstractVector, M<:AbstractMatrix}
    A = convert(M, IdMat(length(s.δ)))
    return Similarity(s.ρ, s.δ, A, s.ρ*A)
end

# Base.promote_rule(Type{Similarity{T, V, M}}, Type{TranslatingSimilarity{T, V}}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix} = Similarity{T, V, M}
Base.promote_rule(::Type{Similarity{R, T, M}}, ::Type{TranslatingSimilarity{R, T}}
                ) where {R<:Real, T<:AbstractVector, M<:AbstractMatrix} = Similarity{R, T, M}

# composition of two similarities
simcomp(s₁::S, s₂::S) where S<:AbstractSimilarity = S(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A*s₂.A, s₁.ρA*s₂.ρA)
simcomp(s₁::TranslatingSimilarity, s₂::Similarity) = Similarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρ*s₂.δ, s₂.A, s₂.ρA)
# (does not commute)
simcomp(s₁::Similarity,s₂::TranslatingSimilarity) = Similarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A, s₁.ρA)
simcomp(s₁::TranslatingSimilarity,s₂::TranslatingSimilarity) = TranslatingSimilarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, one(typeof(s₁.ρ*s₂.ρ)), s₁.ρ*s₂.ρ)
# one dimensional case
# simcomp(s₁::OneDimensionalSimilarity, s₂::OneDimensionalSimilarity) = OneDimensionalSimilarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A*s₂.A, s₁.ρA*s₂.ρA)

# Similarities as maps, optimised as far as possible

# similarity acting on an IFS
function simcompifs(s::T, ifs::AbstractVector{T}) where T<:AbstractSimilarity
    Achain = [s.A*ifs[m].A/s.A for m in eachindex(ifs)]
    rAchain = [ifs[m].ρ*Achain[m] for m in eachindex(ifs)]
    return [T(ifs[m].ρ, (IdMat-rAchain[m])*s.δ + s.ρA*ifs[m].δ, Achain[m], rAchain[m]) for m in eachindex(ifs)]
end

# simplified case without rotation
function simcompifs(s::TranslatingSimilarity, IFS::Vector{<:TranslatingSimilarity})
    return [TranslatingSimilarity(IFS[m].ρ, (IdMat-IFS[m].ρA)*s.δ + s.ρA*IFS[m].δ, IFS[m].A, IFS[m].ρA) for m in eachindex(IFS)]
end

(s::AbstractSimilarity)(x)  = s.ρA*x + s.δ
(s::TranslatingSimilarity)(x) = s.ρ*convert(typeof(s.δ),x) + s.δ

# nested composition:
function simmulticomp(IFS::Vector{<:AbstractSimilarity}, m::AbstractVector{<:Integer})
    s = IFS[m[1]]
        for j in 2:length(m)
            s = simcomp(s,IFS[m[j]])
        end
    return s
end

#inverse map
sim_map_inv(s::Similarity, x) = s.rA/(x-s.δ)

# determine how Similarity appears in the REPL

info_string(s::TranslatingSimilarity) = "x ↦ "*string(round.(s.ρ,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::Similarity) = "x ↦ "*string(round.(s.ρ,digits=2))*string(round.(s.A,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::OneDimensionalSimilarity) = "x ↦ "*string(round(s.ρA,digits=2))*"x + "*string(round(s.δ,digits=2))

Base.show(io::IO, s::AbstractSimilarity)  = print(io,'\n',typeof(s),':','\n', info_string(s))

# fixed points

# x = ρAx + δ
# (I - ρA)x = δ
# x = δ \ (I-ρA)
fixed_point(s::AbstractSimilarity) =  (IdMat - s.ρA) \ s.δ
