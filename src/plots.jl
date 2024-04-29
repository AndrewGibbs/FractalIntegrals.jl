# seems that markers aren't visible if smaller than this:
MIN_MARKER_SIZE = 0.1

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

function plothelper(Γ::AbstractAttractor,
                    max_num_pts,
                    markersize,
                    markersize_const,
                    marker_adjust_scale)

    # get points and convert to x,y coords for plotting
    pts = sketch_attractor(Γ; max_num_pts = max_num_pts)
    xpts = [pt[1] for pt in pts]
    ypts = [pt[2] for pt in pts]

    # default marker size is based on def'n of Hausdorff dimension
    if markersize <= 0
        markersize = marker_adjust_scale*
                    max((markersize_const / length(pts))^(1/Γ.d) *
                            diam(Γ),
                        MIN_MARKER_SIZE)
    end

    return xpts, ypts, markersize
end

function plothelper(μ::AbstractInvariantMeasure,
                max_num_pts,
                markersize,
                markersize_const,
                marker_adjust_scale)

    # get points and convert to x,y coords for plotting
    pts, weights = sketch_measure(μ; max_num_pts = max_num_pts)
    xpts = [pt[1] for pt in pts]
    ypts = [pt[2] for pt in pts]

    # default marker size is based on def'n of Hausdorff dimension
    if markersize <= 0
    markersize = marker_adjust_scale*
        max((markersize_const / length(pts))^(1/μ.supp.d) *
                diam(μ),
            MIN_MARKER_SIZE)
    end

    return xpts, ypts, markersize, weights
end

function Plots.plot( Γ::AbstractAttractor,
                args...;
                max_num_pts::Integer = Integer(1e4),
                markershape = :circle,
                linewidth::Real = 0.0,
                markersize::Real = 0.0,
                markersize_const::Real = 1e3,
                marker_adjust_scale::Real = 1.0,
                label = "",
                kwargs...)
        
        # start by getting main ingredients from helper function:
        xpts, ypts, markersize =
            plothelper( Γ,
                        max_num_pts,
                        markersize,
                        markersize_const,
                        marker_adjust_scale)

        Plots.plot(xpts, ypts;
                markershape = markershape,
                linewidth = linewidth,
                markersize = markersize,
                label = label,
                markerstrokewidth = 0, # no border required
                kwargs...)
end

function Plots.plot!( Γ::AbstractAttractor,
                args...;
                max_num_pts::Integer = Integer(1e5),
                markershape = :circle,
                linewidth::Real = 0.0,
                markersize::Real = 0.0,
                markersize_const::Real = 1e3,
                marker_adjust_scale::Real = 1.0,
                label = "",
                kwargs...)

    # start by getting main ingredients from helper function:
    xpts, ypts, markersize =
        plothelper( Γ,
                    max_num_pts,
                    markersize,
                    markersize_const,
                    marker_adjust_scale)

    Plots.plot!(xpts, ypts;
                markershape = markershape,
                linewidth = linewidth,
                markersize = markersize,
                label = label,
                markerstrokewidth = 0, # no border required
                kwargs...)
end

function Plots.plot( μ::AbstractInvariantMeasure,
                    args...;
                    max_num_pts::Integer = Integer(1e4),
                    markershape = :circle,
                    linewidth::Real = 0.0,
                    markersize::Real = 0.0,
                    markersize_const::Real = 1e3,
                    marker_adjust_scale::Real = 1.0,
                    label = "",
                    kwargs...)

    # start by getting main ingredients from helper function:
    xpts, ypts, markersize, weights =
    plothelper( μ,
                max_num_pts,
                markersize,
                markersize_const,
                marker_adjust_scale)

    Plots.plot(xpts, ypts;
        markershape = markershape,
        linewidth = linewidth,
        markersize = markersize,
        label = label,
        marker_z = weights,
        markerstrokewidth = 0, # no border required
        kwargs...)
end

function Plots.plot!( μ::AbstractInvariantMeasure,
                    args...;
                    max_num_pts::Integer = Integer(1e5),
                    markershape = :circle,
                    linewidth::Real = 0.0,
                    markersize::Real = 0.0,
                    markersize_const::Real = 1e3,
                    marker_adjust_scale::Real = 1.0,
                    label = "",
                    kwargs...)

    # start by getting main ingredients from helper function:
    xpts, ypts, markersize, weights =
    plothelper( μ,
            max_num_pts,
            markersize,
            markersize_const,
            marker_adjust_scale)

    Plots.plot!(xpts, ypts;
        markershape = markershape,
        linewidth = linewidth,
        markersize = markersize,
        label = label,
        marker_z = weights,
        markerstrokewidth = 0, # no border required
        kwargs...)
end