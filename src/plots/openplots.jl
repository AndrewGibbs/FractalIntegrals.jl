# ----------------------- methods for plotting d=2 attractors ----------------- #

# convert my Polygon struct to Shape struct which is used in Plot function
Plots.Shape(p::Polygon) = Shape([x[1] for x in p.nodes], [x[2] for x in p.nodes])

# main plotting function which plots a basis, and can be used to colour the elements
function Plots.plot(basis::QuasiUniformBasis{LebesgueMeasure{2}};
            levels = None,
            colour_map = :jet
            kwargs...)

    # get supporting attractor
    Γ = basis.measure.supp

    # choose number of levels to go down in parent polygon approx
    if levels == None
        default_refinement = 5
        min_depth = minimum.([length(el.vindex) for el in basis])
        levels = max(default_refinement/min_depth, 1)
    end

    # approximate the boundary
    Γ_boundary_approx = polygon_approx(Γ; levels = levels)

    # initialise empty vector of Shapes, which we will define via mesh
    plot_els = Vector{typeof(Shape(Γ_boundary_approx))}(undef, length(basis))

    for (n, el) in enumerate(basis)
        # initialise mesh element 
        pel_polygon = copy(Γ_boundary_approx)

        # apply maps of mesh element in reverse order to get approximation to mesh element bdry
        for mᵢ ∈ reverse(el.vindex)
            pel_polygon = Γ.ifs[mᵢ](pel_polygon)
        end

        # convert polygon to shape type
        plot_els[n] = Plots.Shape(pel_polygon)
    end

    Plots.plot(plot_els;
        c=colour_map, mc=colour_map, fill_z=permutedims(vals),
        labels=:none, kwargs...)
end

# default case warning for use below
function default_polygon_plot_vals(density::Projection)
    @warn "No colouring values provided, using real part of projection coefficients"
    return real(density.coeffs).*[el.normalisation for el in density.basis]
end

# common use case for plotting density
Plots.plot(density::Projection{QuasiUniformBasis{LebesgueMeasure{2}}},
    vals::AbstractVector{<:Real} = default_polygon_plot_vals(density);
    linewidth = 0.0,
    kwargs...
    ) =
    Plots.plot(density.basis;
        fill_z = permutedims(vals),
        linewidth = linewidth,
        kwargs...)