# module FractalPresets
# using .FractalIntegrals

function cantorset(; ρ = 1/3)
    ifs = [Similarity(ρ, 0.0), Similarity(ρ, 1-ρ)]
    reflectiongroup1d = [IdentityInvariantMap(1), OneDimensionalInvariantMap(1, -1)]
    return Attractor(ifs, d = log(2)/log(1/ρ), diam = 1, symmetries = reflectiongroup1d)
end

function cantordust(; ρ = 1/3)
    ifs = [ Similarity(ρ, [0.0, 0.0]),
            Similarity(ρ, [1-ρ, 0.0]),
            Similarity(ρ, [0.0, 1-ρ]),
            Similarity(ρ, [1-ρ, 1-ρ])]
    return Attractor(ifs, d = log(4)/log(1/ρ), diam = sqrt(2), symmetries = DihedralGroup(4))
end

# end