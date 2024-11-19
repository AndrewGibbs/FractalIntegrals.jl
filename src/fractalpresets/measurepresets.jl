function getmeasure(T::Type, fractalname::Union{Symbol, String};
                    suppmeasure = Nothing, kwargs...)

    knownmeasures = Dict(:cantorset => T(1),
                        :kochsnowflake => T(2)*T(sqrt(3))/T(5),
                        :heighwaydragon => T(1/2))
                        
    Γ = getfractal(fractalname; kwargs)

    if isa(Γ, OpenAttractorUnion)
        measuretype = "Lebesgue"
    else
        measuretype = "Hausdorff"
    end

    # get default 
    if suppmeasure == Nothing
        if haskey(knownmeasures, fractalname)
            suppmeasure = knownmeasures(fractalname)
        else
            fractalnamestring = String(fractalname)
            @warn "The $measuretype measure of $fractalnamestring is not known. 
                Thus, the measure will be scaled by the Hausdorff measure of the attractor. 
                Note that this will lead to meaningless results for second-kind integral equations."
        end
    end

    if measuretype == "Lebesgue"
        return LebesgueMeasure(Γ)
    else
        return HausdorffMeasure(Γ)
    end
end

# Option to pass String instead of Symbol. Converts to lowercase and removes spaces from string.
getmeasure(T::Type, fractalname::String; vargs...) =
getmeasure(T::Type, Symbol(lowercase(replace(fractalname, " " => ""))); vargs...)

# define default type
getmeasure(fractalname; vargs...) = getmeasure(Float64, fractalname; vargs...)