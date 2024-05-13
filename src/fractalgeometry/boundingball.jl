function get_boundingball_centre(ifs::AbstractArray{<:AbstractSimilarity})
    M = length(ifs)
    ndims = length(ifs[1].δ)
    𝐀 = sum(ifs[m].ρA for m in 1:M)/M *IdMat(ndims)
    𝐁 = sum(ifs[m].δ for m in 1:M)/M
    # ball centre is given by Appendix A:
    return (IdMat(ndims) - 𝐀) \ 𝐁
end

get_boundingball_centre(Γ::AbstractAttractor) = get_boundingball_centre(Γ.ifs)
get_boundingball_centre(μ::AbstractInvariantMeasure) = get_boundingball_centre(μ.supp.ifs)

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