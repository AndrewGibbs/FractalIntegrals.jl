
function sketch_measure(μ::AbstractInvariantMeasure; max_num_pts::Integer = Integer(1e5))
    submeasures =
        [μ[m] for m in subdivide_indices(μ.supp, 0; max_num_indices = max_num_pts)]

    # return points and measures of subcomponents
    return [get_barycentre(γ) for γ in submeasures], [γ.suppmeasure for γ in submeasures]
end

function sketch_attractor(Γ::AbstractAttractor; max_num_pts::Integer = Integer(1e5))
    return [
        get_boundingball_centre(Γ[𝐦]) for
        𝐦 in subdivide_indices(Γ, 0; max_num_indices = max_num_pts)
    ]
end
