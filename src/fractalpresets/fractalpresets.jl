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
                    # Lebesgue measurable d=2 attractors / dragons
                    :heighwaydragon => heighwaydragon,
                    :kochsnowflake => kochsnowflake,
                    :kochflake => kochsnowflake, # alt name for kochsnowflake
                    :kochcurve => kochcurve,
                    :sierpinskicarpet => sierpinskicarpet,
                    :carpet => sierpinskicarpet # alt name for sierpinskicarpet
                    )

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