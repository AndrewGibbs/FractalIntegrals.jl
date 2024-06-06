# seems that markers aren't visible if smaller than this:
MIN_MARKER_SIZE = 0.1
MAX_PLOT_PTS =  Integer(1e3)

include("sketch.jl")

function plothelper(Γ::AbstractAttractor,
                    max_num_pts,
                    markersize,
                    markersize_const,
                    marker_adjust_scale)

    @assert Γ.n in [1,2] "ambient dimension of attractor must be at most two."
    # get points and convert to x,y coords for plotting
    pts = sketch_attractor(Γ; max_num_pts = max_num_pts)
    xpts = [pt[1] for pt in pts]
    Γ.n == 2 ? ypts = [pt[2] for pt in pts] : ypts = zeros(length(xpts))

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


    @assert μ.supp.n in [1,2] "ambient dimension of support must be at most two."

    # get points and convert to x,y coords for plotting
    pts, weights = sketch_measure(μ; max_num_pts = max_num_pts)
    xpts = [pt[1] for pt in pts]
    μ.supp.n == 2 ? ypts = [pt[2] for pt in pts] : ypts = zeros(length(xpts))

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
                max_num_pts::Integer = MAX_PLOT_PTS,
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
                max_num_pts::Integer = MAX_PLOT_PTS,
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
                    max_num_pts::Integer = MAX_PLOT_PTS,
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
                    max_num_pts::Integer = MAX_PLOT_PTS,
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