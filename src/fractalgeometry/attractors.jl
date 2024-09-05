# ------------------ definitions of different attractor types -----------------------#

abstract type AbstractAttractor{N, M, T} end
abstract type AbstractHomogenousAttractor{N, M, T} <: AbstractAttractor{N, M, T} end

struct Attractor{N, M, T} <: AbstractAttractor{N, M, T}
    ifs::SVector{M, Similarity{N, T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N, T}}
end

struct HomogenousAttractor{N, M, T} <: AbstractHomogenousAttractor{N, M, T}
    ifs::SVector{M, Similarity{N, M, T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{Similarity{N, T}}
    ρ::T
end

struct OneDimensionalAttractor{M, T} <: AbstractAttractor{1, M, T}
    ifs::SVector{M, OneDimensionalSimilarity{N, T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{OneDimensionalSimilarity{T}}
end

struct OneDimensionalHomogenousAttractor{M, T} <: AbstractHomogenousAttractor{1, M, T}
    ifs::SVector{M, OneDimensionalSimilarity{N, T}}
    diam::T
    d::T
    connectedness::Matrix{Bool}
    symmetries::Vector{OneDimensionalSimilarity{T}}
    ρ::T
end

# struct HomogenousNonRotatingAttractor{  T,
#                                         R<:Real,
#                                         Sifs<:AbstractArray{<:AbstractSimilarity{R, T}},
#                                         Ssym<:AbstractArray{<:AbstractSimilarity{R, T}}
#                                         } <: AbstractHomogenousAttractor{R, T}
#     ifs::Sifs
#     diam::R
#     d::R
#     n::Int64
#     connectedness::Matrix{Bool}
#     symmetries::Ssym
#     ρ::R
# end


# define eltype for attractors - will be useful elsewhere
Base.eltype(::AbstractAttractor{N, M, T}) where {N, M, T} = T

# ----------- mappings of the underlying IFS to vector points --------------- 

function ifs_map!(  Sx::AbstractVector{<:T},
                    S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{<:T}) where T
    # M = length(S)
    N = length(x)
    # y = zeros(T, M*N)
    for m in eachindex(S)
        Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
    end
end

function ifs_map(   S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{T}) where T
    N = length(x)
    Sx = zeros(T, length(S)*N)
    for m in eachindex(S)
        @inbounds Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
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

    for 𝐦 = index[end:-1:1]
        new_diam *= Γ.ifs[𝐦].ρ
        new_ifs = simcompifs(Γ.ifs[𝐦], new_ifs)
        new_symmetries = simcompifs(Γ.ifs[𝐦], new_symmetries)
        # new_symmetries = simcompsymmetries(Γ.ifs[m], new_symmetries)
    end

    return new_diam, new_ifs, new_symmetries
end

function get_subattractor(Γ::HomogenousAttractor, 𝐦::AbstractVector{<:Integer})
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, 𝐦)
    return HomogenousAttractor(new_ifs,
            new_diam,
            Γ.d,
            Γ.connectedness,
            new_symmetries,
            Γ.ρ)
end

function get_subattractor(Γ::A, 𝐦::AbstractVector{<:Integer}) where A<:AbstractAttractor
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, 𝐦)
    return A(new_ifs,
            new_diam,
            Γ.d,
            Γ.n,
            Γ.connectedness,
            new_symmetries)
end

# -------------------- outer constructor ---------------------------------- #
function ishomogeneous(ifs::AbstractVector{<:AbstractSimilarity})
    # check if fractal is homogeneous
    homogeneous = true
    for j in 1:(length(ifs)-1)
        if  !(ifs[j].ρ ≈ ifs[j+1].ρ)
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
        d = log(1/length(r))/log(r[1])
    else
        # approximate Hausdorff dimension by approximating zero of following:
        f(d) = sum(r.^d) - 1
        # over range (0,n]
        d = find_zero(f, (0, n*(1+10*eps(ifs[1].ρ))), Bisection())
    end
    return d
end

# user-friendly constructor
function Attractor( ifs::AbstractVector{S};
                    diam = diam(ifs),
                    d::Real = dimH(ifs), # integers welcome
                    connectedness = IdMat(N),
                    symmetries = trivialgroup(T, N)
                    ) where {N, T, S<:AbstractSimilarity{N, T}}

                    
    diamT, dT, _ = promote(diam, d, T) # ensure these are of same type
    # Not type-stable. But I don't think this will be a performance-critical function in practice.
    if N == 1
        if ishomogeneous(ifs)
            Γ  = OneDimensionalHomogenousAttractor(
                    ifs,
                    diamT,
                    dT,
                    connectedness,
                    symmetries,
                    ifs[1].ρ
                    )
                    
        else
            Γ  = OneDimensionalAttractor(
                    ifs,
                    diamT,
                    dT,
                    connectedness,
                    symmetries
                    )
        end
    else
        if ishomogeneous(ifs)
            Γ  = HomogenousAttractor(
                    ifs,
                    diamT,
                    dT,
                    connectedness,
                    symmetries,
                    ifs[1].ρ
                    )
                    
        else
            Γ  = Attractor(
                    ifs,
                    diamT,
                    dT,
                    connectedness,
                    symmetries
                    )
        end
    end
    return Γ
end