# ------------------------ functions for bounding ball ------------------------

function get_boundingball_centre(ifs::AbstractArray{<:AbstractSimilarity})
    M = length(ifs)
    ndims = length(ifs[1].δ)
    𝐀 = sum(ifs[m].ρA for m in 1:M)/M *IdMat(ndims)
    𝐁 = sum(ifs[m].δ for m in 1:M)/M
    return (IdMat(ndims) - 𝐀) \ 𝐁
end

function get_boundingball_centre(ifs::AbstractArray{OneDimensionalSimilarity{T}}
                                ) where {T<:Real}
    min_pt = T(Inf)
    max_pt = -T(Inf)
    # consider extremal fixed points
    for s in ifs
        min_pt = min(fixed_point(s),min_pt)
        max_pt = max(fixed_point(s),max_pt)
    end
    
    # take halfway point between endpoints
    return (min_pt + max_pt)/2
end

get_boundingball_centre(Γ::AbstractAttractor) = get_boundingball_centre(Γ.ifs)
get_boundingball_centre(μ::AbstractInvariantMeasure) = get_boundingball_centre(μ.supp.ifs)

# ---------------------- diameter of attractor given IFS ------------------------------
function diam(ifs::AbstractArray{<:AbstractSimilarity})
    M = length(ifs)
    # uses idea of `Spatial Bounding of Self-Affine Iterated Function System Attractor Sets'
    # Jonathon Rice, Graphics interface, 1996. Computes bounding ball, diameter is inferred.
    # ball centre is given by Appendix A:
    ball_centre = get_boundingball_centre(ifs)
    # diameter is given by eqn (12):
    radius = maximum(norm(ball_centre - ifs[m](ball_centre))/(1-ifs[m].ρ) for m in 1:M)
    return 2*radius
end

diam(Γ::AbstractAttractor) = Γ.diam
diam(μ::AbstractInvariantMeasure) = μ.supp.diam

# ---------------------- distance from an attractor to a point ---;; ----------------#
dist⁺(Γ::AbstractAttractor, x) = sqrt(norm(x - get_boundingball_centre(Γ))^2 + (diam(Γ)/2)^2)
dist⁻(Γ::AbstractAttractor, x) = max(norm(x - get_boundingball_centre(Γ)) - diam(Γ)/2, 0.0)

# -------------------- distance from one attractor to another ------------------ #
dist⁺(Γ::AbstractAttractor, γ::AbstractAttractor) = norm(get_boundingball_centre(γ) - get_boundingball_centre(Γ))
dist⁻(Γ::AbstractAttractor, γ::AbstractAttractor) = dist⁺(Γ, γ) - diam(Γ)/2 - diam(γ)/2