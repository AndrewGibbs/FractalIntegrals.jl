abstract type AbstractAttractor{T, R<:Real} end
abstract type AbstractHomogenousAttractor{T, R} <: AbstractAttractor{T, R} end

struct HomogenousAttractor{ T,
                            R<:Real,
                            S<:AbstractArray{<:AbstractSimilarity{R, T}},
                            G<:AbstractArray{<:AbstractInvariantMap{T}}
                            } <: AbstractHomogenousAttractor{T, R}
    ifs::S
    diam::R
    d::R
    n::Int64
    connectedness::Matrix{Bool}
    symmetries::G
    ρ::R
end

struct HomogenousNonRotatingAttractor{  T,
                                        R<:Real,
                                        S<:AbstractArray{TranslatingSimilarity{R, T}},
                                        G<:AbstractArray{<:AbstractInvariantMap{T}}
                                        } <: AbstractHomogenousAttractor{T, R}
    ifs::S
    diam::R
    d::R
    n::Int64
    connectedness::Matrix{Bool}
    symmetries::G
    ρ::R
end

struct Attractor{   T,
                    R<:Real,
                    S<:AbstractArray{<:AbstractSimilarity{R, T}},
                    G<:AbstractArray{<:AbstractInvariantMap{T}}
                    } <: AbstractAttractor{T, R}
    ifs::S
    diam::R
    d::R
    n::Int64
    connectedness::Matrix{Bool}
    symmetries::G
end

# ALSO DEFINE HomogenousNonRotatingAttractor, OneDimensionalAttractor and OneDimensionalHomogenousAttractor
# - this will be very useful and elegant later on
# the most inelegand thing about it is there is no type hierarcy - it's a type matix :/
# but this can be addressed by taking type unions in method calls,
# and I can choose a heirarcy that fits best - I think AbstractHomogenousAttractor < AbstractAttractor?

# Actually based on this:
#   AbstractInvariantMeasure{<:AbstractAttractor{T, R}} where {T <: Real, R}
# things will still be quite elegant without the need for OneDimensionalBlahBlah
# which means I can have a heirarchy easily

# define eltype for attractors - will be useful elsewhere
Base.eltype(::AbstractAttractor{T, R}) where {T, R} = T

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

function get_subattractor_elements(Γ::AbstractAttractor, index::AbstractVector{<:Integer})
    # get new measure and diameter. First initialise:
    new_diam = Γ.diam
    new_ifs = Γ.ifs
    new_symmetries = Γ.symmetries

    for m = index[end:-1:1]
        new_diam *= Γ.ifs[m].ρ
        new_ifs = simcompifs(Γ.ifs[m], new_ifs)
        new_symmetries = simcompsymmetries(Γ.ifs[m], new_symmetries)
    end

    return new_diam, new_ifs, new_symmetries
end

function get_subattractor(Γ::HomogenousAttractor, m::AbstractVector{<:Integer})
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, m)
    return HomogenousAttractor(new_ifs,
            new_diam,
            Γ.d,
            Γ.n,
            Γ.connectedness,
            new_symmetries,
            Γ.ρ)
end

function get_subattractor(Γ::A, m::AbstractVector{<:Integer}) where A<:Attractor
    new_diam, new_ifs, new_symmetries = get_subattractor_elements(Γ, m)
    return A(new_ifs,
            new_diam,
            Γ.d,
            Γ.n,
            Γ.connectedness,
            new_symmetries)
end

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
                    connectedness = Matrix(IdMat(length(ifs))),
                    symmetries = trivialgroup(typeof(ifs[1].ρ), length(ifs[1].δ))
                    ) where S<:AbstractSimilarity

    n = length(ifs[1].δ) # ambient dimensions
    diamT, dT = promote(diam, d) # ensure these are of same type
    # Not type-stable. But I don't think this will be a performance-critical function in practice.
    if isa(S, AbstractSimilarity)
        Γ  = HomogenousNonRotatingAttractor(
                ifs,
                diamT,
                dT,
                n,
                connectedness,
                symmetries,
                ifs[1].ρ
                )

    elseif ishomogeneous(ifs)
        Γ  = HomogenousAttractor(
                ifs,
                diamT,
                dT,
                n,
                connectedness,
                symmetries,
                ifs[1].ρ
                )
                
    else
        Γ  = Attractor(
                ifs,
                diamT,
                dT,
                n,
                connectedness,
                symmetries
                )
    end
    return Γ
end