abstract type AbstractSimilarity{N, T<:Number} end
# abstract type VectorSimilarity{N,T} <: AbstractSimilarity{T} end

# ----------------- concrete similarity types -----------------------#

"""
    Similarity{N, T} <: AbstractSimilarity{N, T}
Constructs a similarity, which can be used to describe geometry of self-similar fractals.
`N` is the dimension of the ambient space, and `T<:Number`.

Mathematically, this captures the following map:

```math
s(x) = δ + ρAx, x ∈ ℝⁿ
```

# Fields
- `ρ :: T`: The contraction factor.
- `δ :: SVector{N,T}`: The translation vector.
- `A :: SMatrix{N,N,T}` : The rotation/reflection matrix.
- `ρA` :: SMatrix{N,N,T} : The rotation/reflection matrix multiplied by the contraction factor.

# Examples
```julia
s₁ = Similarity(1/3, [0,0])
s₂ = Similarity(Float32(1/3), rand(3), [-1 0; 0 1])
```
In the first example above, the third argument is not given, 
so it is assumed that `A=I`.

Various binary operations can be performed on Similarity:
```julia
s₂₁ = s₁ ∘ s₂
s₂⁻¹ = inv(s₂)
s₁³ = s₁^3
```
In each of the above examples, a Similarity is returned.

Similarity can be applied as a function
```julia
x = rand(2)
y = s₁(x)
```
"""
struct Similarity{N, T} <: AbstractSimilarity{N, T}# <: VectorSimilarity{N,T}
    ρ :: T
    δ :: SVector{N,T}
    A :: SMatrix{N,N,T}
    ρA :: SMatrix{N,N,T}
end

# struct TranslatingSimilarity{N,T} <: VectorSimilarity{N,T}
#     ρ :: T
#     δ :: SVector{N,T}
#     A :: T
#     ρA :: T
# end

struct OneDimensionalSimilarity{T} <: AbstractSimilarity{1, T}
    ρ :: T
    δ :: T
    A :: T
    ρA :: T
end

# ---------------- outer constructors for similarity types  ---------------- #

# function Similarity(ρ::T1, δ::AbstractVector{T2}) where {T1<:Number, T2<:Number}
#     # @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
#     T = promote_type(T1, T2)
#     return TranslatingSimilarity(T(ρ), SVector{length(δ),T}(δ), one(T), T(ρ))
# end

function Similarity(ρ::T1, δ::AbstractVector{T2}, A::AbstractMatrix{T3}=IdMat(length(δ))
                    ) where {T1<:Number, T2<:Number, T3<:Number}
    # @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(det(A)) ≈ 1 "Matrix (third argument) must have determinant one"
    N = length(δ)
    @assert (N,N) == size(A) "Matrix dimension must match translation vector"
    T = promote_type(T1,T2,T3)
    return Similarity(T(ρ), SVector{N,T}(δ), SMatrix{N,N,T}(A), SMatrix{N,N,T}(ρ*A))
end

function Similarity(ρ::T1, δ::T2, r::T3=one(T1)) where {T1<:Number, T2<:Number, T3<:Number}
    # @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(r) ≈ 1 "Reflection (third argument) must have length one"
    T = promote_type(T1,T2,T3)
    return OneDimensionalSimilarity(T(ρ), T(δ), T(r), convert(T,ρ*r))
end

# ------------------- composition of similarities -------------------------------------- #

simcomp(s₁::S, s₂::S) where S<:AbstractSimilarity = S(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A*s₂.A, s₁.ρA*s₂.ρA)

# overload composition operator for syntactic sugar
Base.:∘(s₁::AbstractSimilarity, s₂::AbstractSimilarity) = simcomp(s₁,s₂)

Base.inv(s::S) where S<:AbstractSimilarity = S(inv(s.ρ), -s.ρA\s.δ, inv(s.A), inv(s.ρA))

function Base.:^(s::AbstractSimilarity, p::Integer)
    if p<0
        s⁻¹ = inv(s)
        sᵖ = (s⁻¹)^(-p) # call recursively (once)
    else
        sᵖ = s
        # apply repeated composition
        for _= 2:p
            sᵖ = sᵖ ∘ s
        end
    end
    return sᵖ
end

# -------------------------------- Similarities as maps --------------------------------------

(s::AbstractSimilarity)(x)  = s.ρA*x + s.δ

# similarity acting on an iterated function system (IFS)
function simcompifs(s::AbstractSimilarity, ifs::AbstractVector{T}) where T<:AbstractSimilarity
    Achain = [s.A*ifs[m].A/s.A for m in eachindex(ifs)]
    rAchain = [ifs[m].ρ*Achain[m] for m in eachindex(ifs)]
    return [T(ifs[m].ρ, (IdMat-rAchain[m])*s.δ + s.ρA*ifs[m].δ, Achain[m], rAchain[m]) for m in eachindex(ifs)]
end

# nested composition:
function simmulticomp(IFS::AbstractVector{<:AbstractSimilarity}, m::AbstractVector{<:Integer})
    s = IFS[m[1]]
        for j in 2:length(m)
            s = simcomp(s,IFS[m[j]])
        end
    return s
end

# ------------------------------------------- inverse map --------------------------#
sim_map_inv(s::Similarity, x) = s.rA/(x-s.δ)

# ------------------ determine how Similarity appears in the REPL ------------------#

info_string(s::Similarity) = "x ↦ "*string(round.(s.ρ,digits=2))*string(round.(s.A,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::OneDimensionalSimilarity) = "x ↦ "*string(round(s.ρA,digits=2))*"x + "*string(round(s.δ,digits=2))

Base.show(io::IO, s::AbstractSimilarity)  = print(io,'\n',typeof(s),':','\n', info_string(s))

# ------------------------------ fixed points --------------------------------------#
fixed_point(s::AbstractSimilarity) =  (IdMat - s.ρA) \ s.δ


# ------------------------ define identity and zero similarities ---------------------------- #

Base.one(::OneDimensionalSimilarity{T}) where {T} =
    OneDimensionalSimilarity{T}(one(T), zero(T), one(T), one(T))

Base.zero(::OneDimensionalSimilarity{T}) where {T} =
    OneDimensionalSimilarity{T}(zero(T), zero(T), zero(T), zero(T))

Base.one(::Similarity{N, T}) where {N, T} =
    Similarity{N, T}(one(T), SVector{N, T}(zeros(T, N)), SMatrix{N, N}(Matrix{T}(IdMat(N))), SMatrix{N, N}(Matrix{T}(IdMat(N))))

Base.zero(::Similarity{N, T}) where {N, T} =
    Similarity{N, T}(zero(T), SVector{N, T}(zeros(T, N)), SMatrix{N, N}(zeros(T, N, N)), SMatrix{N, N}(zeros(T, N, N)))