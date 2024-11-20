# ----------------------- methods for plotting d=2 attractors ----------------- #

# convert my Polygon struct to Shape struct which is used in Plot function
Plots.Shape(p::Polygon) = Shape([x[1] for x in p.nodes], [x[2] for x in p.nodes])

# main plotting function which plots a basis, and can be used to colour the elements
function Plots.plot(basis::QuasiUniformBasis{<:LebesgueMeasure{2}},
                    vals::AbstractVector{<:Real};
                    default_refinement = 5,
                    levels = Nothing,
                    colour_map = :jet,
                    kwargs...)

    # get supporting attractor
    Γ = basis.measure.supp

    # choose number of levels to go down in parent polygon approx
    if levels == Nothing
        min_depth = minimum([length(el.vindex) for el in basis])
        levels = max(default_refinement-min_depth, 1)
    end

    # approximate the boundary
    Γ_boundary_approx = polygon_approx(Γ; levels = levels)

    # initialise empty vector of Shapes, which we will define via mesh
    plot_els = Vector{typeof(Shape(Γ_boundary_approx))}(undef, length(basis))

    for (n, el) in enumerate(basis)
        # # initialise mesh element 
        pel_polygon = Γ_boundary_approx

        # apply maps of mesh element in reverse order to get approximation to mesh element bdry
        for mᵢ ∈ reverse(el.vindex)
            pel_polygon = Γ.ifs[mᵢ](pel_polygon)
        end

        # convert polygon to shape type
        plot_els[n] = Plots.Shape(pel_polygon)
    end

    Plots.plot(plot_els;
        fill_z = permutedims(vals),
        c=colour_map, mc=colour_map,
        labels=:none, kwargs...)
end

# default case warning for use below
function get_polygon_proj_vals(density::Projection, plot_type)
    denormalisation = [el.normalisation for el in density.basis]
    if plot_type == :none
        @warn "No colouring values provided, using real part of projection coefficients"
        return real(density.coeffs).*denormalisation
    elseif plot_type == :real
        return real(density.coeffs).*denormalisation
    elseif plot_type == :imag
        return imag(density.coeffs).*denormalisation
    elseif plot_type == :abs
        return abs.(density.coeffs).*denormalisation
    end
end

# common use case for plotting density
function Plots.plot(density::Projection{<:QuasiUniformBasis{<:LebesgueMeasure{2}}};
    plot_type = :none,
    linewidth = 0.0,
    kwargs...
    )
    vals = get_polygon_proj_vals(density, plot_type)
    Plots.plot( density.basis,
                vals;
                linewidth = linewidth,
                kwargs...)

end