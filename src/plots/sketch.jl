
function sketch_measure(μ::AbstractInvariantMeasure; max_num_pts::Integer = Integer(1e5))

    submeasures = [μ[m] for m in subdivide_indices(μ.supp, 0, max_num_pts)]

    # return points and measures of subcomponents
    return [γ.barycentre for γ in submeasures], [γ.suppmeasure for γ in submeasures]
end

function sketch_attractor(Γ::AbstractAttractor; max_num_pts::Integer = Integer(1e5))

    # convert to Hausdorff measure so we can easily get a notion of 'centre' point
    μ = HausdorffMeasure(Γ)

    # get submeasures, equivalent to subdividing attractor
    submeasures = [μ[m] for m in subdivide_indices(μ.supp, 0, max_num_pts)]

    # return barycentres
    return [γ.barycentre for γ in submeasures]
end