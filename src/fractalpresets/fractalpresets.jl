include("compactattractors.jl")
include("dragonsandsnowflakes.jl")

# ---------------------------------------------------------------------------------------- #
# --------------------- main functions to export fractals  ------------------------------- #
# ---------------------------------------------------------------------------------------- #
# the following dict will be used to get fractals without exporting all of these functions
fractaldict = Dict( 
                    :cantorset => cantorset,
                    # measure zero Γ ⊂ ℝ² attractors
                    :cantordust => cantordust,
                    :sierpinskitriangle => sierpinskitriangle,
                    :sierpinski => sierpinskitriangle, # alt name for sierpinskitriangle
                    :sierpinskigasket => sierpinskitriangle, # alt name for sierpinskitriangle
                    :vicsek => vicsek,
                    :kochcurve => kochcurve,
                    :sierpinskicarpet => sierpinskicarpet,
                    :carpet => sierpinskicarpet, # alt name for sierpinskicarpet
                    # Lebesgue measurable d=2 attractors / dragons
                    :heighwaydragon => heighwaydragon,
                    :kochsnowflake => kochsnowflake,
                    :kochflake => kochsnowflake # alt name for kochsnowflake
                    )

"""
    getfractal(T::Type, fractalname::Symbol; vargs...)
    getfractal(T::Type, fractalname::String; vargs...)
    getfractal(fractalname; vargs...)

Returns a preset fractal attractor corresponding to `fractalname`.
This can be a sybol or a string.

# Presets:
- `cantorset`: Cantor Set
- `cantordust`: Cantor Dust
- `sierpinskitriangle`: Sierpinski Triangle
- `vicsek`: Vicsek Fractal
- `kochcurve`: Koch Curve
- `sierpinskicarpet`: Sierpinski Carpet
- `heighwaydragon`: Heighway Dragon
- `kochsnowflake`: Koch Snowflake

Some of the presets will take an optional input argument.
For example, the contraction factor ρ∈(0,1/2] of Cantor sets and dust.
```julia
Γ = getfractal("cantor set"; ρ = 1/4)
```
"""
function getfractal(T::Type, fractalname::Symbol; vargs...)
    if haskey(fractaldict, fractalname)
        return fractaldict[fractalname](T; vargs...)
    else
        str_fractalname = String(fractalname)
        error("Invalid fractal name: $str_fractalname")
    end
end

# Option to pass String instead of Symbol. Converts to lowercase and removes spaces from string.
getfractal(T::Type, fractalname::String; vargs...) =
    getfractal(T::Type, Symbol(lowercase(replace(fractalname, " " => ""))); vargs...)

# define default type
getfractal(fractalname; vargs...) = getfractal(Float64, fractalname; vargs...)


