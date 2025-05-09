# ------------------ definitions of different attractor types -----------------------#
abstract type Fractal end
abstract type AbstractAttractor{N,M,T} <: Fractal end
abstract type AbstractHomogenousAttractor{N,M,T} <: AbstractAttractor{N,M,T} end

"""
    Attractor{N, M, T} <: AbstractAttractor{N, M, T}
An attractor of an iterated function system (IFS).
In this context, an IFS should be interpreted as a set of similarities.
Mathematically, an IFS attractor is the unique bounded non-empty set satisfying
```math
Γ = ⋃ₘ sₘ(Γ)
```
where sₘ are similarities, defined in the docstring `Similarity`.
For a fuller explanation, see for e.g.
this wikipedia article `https://en.wikipedia.org/wiki/Iterated_function_system`.

# Parameters
- `N` is the ambient dimension of the attractor
- `M` is the number of similarities
- `T<:Number` is the numeric type

# Fields
- `ifs::SVector{M, Similarity{N, T}}`: iterated function system
- `diam::T`: diamter of attractor
- `d::T`: Hausdorff dimension of attractor
- `connectedness::Matrix{Bool}`: matrix describing connected subcomponents
- `symmetries::Vector{Similarity{N, T}}`: maps in symmetry group

Attractors can be constructed using a vector or tuple of similarities to represent the IFS;
the diameter and Hausdorff dimension are computed automatically, 
the connectedness matrix is taken to be the identity, corresponding to a disjoint attractor,
and the symmetry group is assumed to be the trivial group.

If the connectedness matrix is assumed to be the identity when the attractor is non-disjoint,
this will impact accuracy of calculations of singular integrals.

# Examples
```julia
s₁ = Similarity(1/3, 0)
s₂ = Similarity(1/3, 2/3)
cantor_set = Attractor(s₁, s₂)

courage = Similarity(1/2, [0, 0])
wisdom  = Similarity(1/2, [1/2, 0])
power   = Similarity(1/2, [1/4, sqrt(3)/4])
sierpinski_triangle = Attractor(courage, wisdom, power, connectedness = Bool(ones(3,3)))
```
Note how the connectedness matrix for the Sierpinski triangle is all ones,
because the subcomponents all touch at the corners.
"""
struct Attractor{N,M,T} <: AbstractAttractor{N,M,T}
    ifs::SVector{M,Similarity{N,T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N,T}}
end

struct OpenAttractor{N,M,T} <: AbstractAttractor{N,M,T}
    ifs::SVector{M,Similarity{N,T}}
    diam::T
    d::Int64
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N,T}}

    # Inner constructor to enforce d == N
    # function OpenAttractor{N, M, T}(ifs::SVector{M, Similarity{N, T}}, diam::T,
    #     d, connectedness::Matrix{Bool}, symmetries::Vector{Similarity{N, T}}
    #     ) where {N, M, T}
    #     if d != N || !isa(d, Integer)
    #         throw(ArgumentError("The field `d` must be an integer equal to parameter `N` (d=$d, N=$N)"))
    #     else
    #         new{N, M, T}(ifs, diam, Int64(d), connectedness, symmetries)
    #     end
    # end
end

struct HomogenousAttractor{N,M,T} <: AbstractHomogenousAttractor{N,M,T}
    ifs::SVector{M,Similarity{N,T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N,T}}
    ρ::T
end

struct OpenHomogenousAttractor{N,M,T} <: AbstractHomogenousAttractor{N,M,T}
    ifs::SVector{M,Similarity{N,T}}
    diam::T
    d::Int64
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N,T}}
    ρ::T

    # Inner constructor to enforce d == N
    # function OpenHomogenousAttractor{N, M, T}(ifs::SVector{M, Similarity{N, T}}, diam::T,
    #     d, connectedness::Matrix{Bool}, symmetries::Vector{Similarity{N, T}}, ρ::T
    #     ) where {N, M, T}
    #     if d != N || !isa(d, Integer)
    #         throw(ArgumentError("The field `d` must be an integer equal to parameter `N` (d=$d, N=$N)"))
    #     else
    #         new{N, M, T}(ifs, diam, Int64(d), connectedness, symmetries, ρ)
    #     end
    # end
end

struct OneDimensionalAttractor{M,T} <: AbstractAttractor{1,M,T}
    ifs::SVector{M,OneDimensionalSimilarity{T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{OneDimensionalSimilarity{T}}
end

struct OneDimensionalHomogenousAttractor{M,T} <: AbstractHomogenousAttractor{1,M,T}
    ifs::SVector{M,OneDimensionalSimilarity{T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{OneDimensionalSimilarity{T}}
    ρ::T
end

# # useful to consider this union later on:
# OneDimensionalAttractorUnion = Union{OneDimensionalAttractor, OneDimensionalHomogenousAttractor}

# # open attractors also need a union
# OpenAttractorUnion = Union{OneDimensionalAttractor, OneDimensionalHomogenousAttractor}

# useful to consider this union later on:
OneDimensionalAttractorUnion{M,T} =
    Union{OneDimensionalAttractor{M,T},OneDimensionalHomogenousAttractor{M,T}}

# open attractors also need a union
OpenAttractorUnion{N,M,T} = Union{OpenAttractor{N,M,T},OpenHomogenousAttractor{N,M,T}}

# define eltype for attractors - will be useful elsewhere
# Base.eltype(::AbstractAttractor{N, M, T}) where {N, M, T} = T
Base.eltype(Γ::AbstractAttractor) = typeof(Γ.ifs[1].δ)

# ----------- mappings of the underlying IFS to vector points --------------- 

function ifs_map!(
    Sx::AbstractVector{<:T},
    S::AbstractVector{<:AbstractSimilarity},
    x::AbstractVector{<:T},
) where {T}
    # M = length(S)
    N = length(x)
    # y = zeros(T, M*N)
    for m in eachindex(S)
        Sx[((m-1)*N+1):(m*N)] .= S[m].(x)
    end
end

function ifs_map(S::AbstractVector{<:AbstractSimilarity}, x::AbstractVector{T}) where {T}
    N = length(x)
    Sx = zeros(T, length(S) * N)
    for m in eachindex(S)
        @inbounds Sx[((m-1)*N+1):(m*N)] .= S[m].(x)
    end
    return Sx
end

(Γ::AbstractAttractor)(x::AbstractVector) = ifs_map(Γ.ifs, x)

# ------------------------ sub-attractors ----------------------------------------#

function get_subattractor_elements(Γ::AbstractAttractor, index::AbstractVector{<:Integer})
    # get new measure and diameter. First initialise:
    new_diam = Γ.diam
    new_ifs = Γ.ifs
    new_symmetries = Γ.symmetries

    for 𝐦 in index[end:-1:1]
        new_diam *= Γ.ifs[𝐦].ρ
        new_ifs = simcompifs(Γ.ifs[𝐦], new_ifs)
        new_symmetries = simcompifs(Γ.ifs[𝐦], new_symmetries)
        # new_symmetries = simcompsymmetries(Γ.ifs[m], new_symmetries)
    end

    return new_diam, new_ifs, new_symmetries
end

function get_subattractor(
    Γ::A,
    𝐦::AbstractVector{<:Integer},
) where {A<:AbstractHomogenousAttractor}
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, 𝐦)
    return A(new_ifs, new_diam, Γ.d, Γ.connectedness, new_symmetries, Γ.ρ)
end

function get_subattractor(Γ::A, 𝐦::AbstractVector{<:Integer}) where {A<:AbstractAttractor}
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, 𝐦)
    return A(new_ifs, new_diam, Γ.d, Γ.connectedness, new_symmetries)
end

# -------------------- outer constructor ---------------------------------- #
function ishomogeneous(ifs::AbstractVector{<:AbstractSimilarity})
    # check if fractal is homogeneous
    homogeneous = true
    for j = 1:(length(ifs)-1)
        if !(ifs[j].ρ ≈ ifs[j+1].ρ)
            homogeneous = false
            break
        end
    end
    return homogeneous
end

function dimH(ifs::AbstractVector{<:AbstractSimilarity})
    n = length(ifs[1].δ)
    r = [s.ρ for s in ifs]

    if ishomogeneous(ifs)
        d = log(1 / length(r)) / log(r[1])
    else
        # approximate Hausdorff dimension by approximating zero of following:
        f(d) = sum(r .^ d) - 1
        # over range (0,n]
        d = find_zero(f, (0, n * (1 + 10 * eps(ifs[1].ρ))), Bisection())
    end
    return d
end

# for completeness
dimH(Γ::AbstractAttractor) = Γ.d

# user-friendly constructor
function Attractor(
    ifs::AbstractVector{S};
    diam = diam(ifs),
    d::Real = dimH(ifs),
    connectedness = Matrix(IdMat(length(ifs))),
    symmetries = trivialgroup(T, N),
) where {N,T,S<:AbstractSimilarity{N,T}}
    diamT, dT = promote(T(diam), T(d)) # ensure these are of same type
    M = length(ifs)
    sv_ifs = SVector{M}(ifs)
    # connectedness = Bool.(connectedness)
    # Not type-stable. But I don't think this will be a performance-critical function in practice.
    if N == 1
        if ishomogeneous(ifs)
            Γ = OneDimensionalHomogenousAttractor(
                sv_ifs,
                diamT,
                dT,
                connectedness,
                symmetries,
                ifs[1].ρ,
            )

        else
            Γ = OneDimensionalAttractor(sv_ifs, diamT, dT, connectedness, symmetries)
        end
    elseif N == d && isa(d, Integer)
        if ishomogeneous(ifs)
            Γ = OpenHomogenousAttractor(
                sv_ifs,
                diamT,
                d,
                connectedness,
                symmetries,
                ifs[1].ρ,
            )

        else
            Γ = OpenAttractor(sv_ifs, diamT, d, connectedness, symmetries)
        end
    else
        if ishomogeneous(ifs)
            Γ = HomogenousAttractor(sv_ifs, diamT, dT, connectedness, symmetries, ifs[1].ρ)

        else
            Γ = Attractor(sv_ifs, diamT, dT, connectedness, symmetries)
        end
    end
    return Γ
end

Attractor(args...; kwargs...) = Attractor([args...]; kwargs...)

# ------------------ determine how Similarity appears in the REPL ------------------#

function Base.show(io::IO, Γ::AbstractAttractor{N,M,T}) where {N,M,T}
    print(io, round(Γ.d; digits = 2), "-dimensional ", typeof(Γ), ":")
    for s in Γ.ifs
        print(io, '\n', info_string(s))
    end
end
