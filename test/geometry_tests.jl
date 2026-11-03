# Properties which must hold for fractals, whatever the implementation

const ALL_PRESETS = [
    "cantorset",
    "cantordust",
    "sierpinskitriangle",
    "sierpinskicarpet",
    "vicsek",
    "kochcurve",
    "kochsnowflake",
    "heighwaydragon",
]

@testset "similarities" begin
    s = Similarity(1 / 3, [0.2, 0.1], [0.0 -1.0; 1.0 0.0])
    t = Similarity(1 / 2, [-0.4, 0.3])
    x = [0.3, 0.7]

    @test s(x) ≈ s.δ + s.ρ * s.A * x
    @test (s ∘ t)(x) ≈ s(t(x))
    @test inv(s)(s(x)) ≈ x
    @test (s^3)(x) ≈ s(s(s(x)))
    @test (s^-1)(s(x)) ≈ x
    @test one(s)(x) ≈ x
    @test norm(s(x) - s([1.0, 2.0])) ≈ s.ρ * norm(x - [1.0, 2.0])

    s₁ = Similarity(1 / 3, 2 / 3)
    @test s₁ isa F.OneDimensionalSimilarity
    @test s₁(1.0) ≈ 1.0
    @test (s₁ ∘ s₁)(0.0) ≈ s₁(s₁(0.0))
    @test inv(s₁)(s₁(0.4)) ≈ 0.4
end

@testset "attractor types from the Attractor constructor" begin
    @test Attractor(Similarity(1 / 3, 0.0), Similarity(1 / 3, 2 / 3)) isa
          F.OneDimensionalHomogenousAttractor
    @test Attractor(Similarity(1 / 3, 0.0), Similarity(1 / 2, 1 / 2)) isa
          F.OneDimensionalAttractor
    @test Attractor(Similarity(1 / 3, [0.0, 0.0]), Similarity(1 / 3, [2 / 3, 0.0])) isa
          F.HomogenousAttractor
    @test Attractor(Similarity(1 / 3, [0.0, 0.0]), Similarity(1 / 2, [1 / 2, 0.0])) isa
          F.Attractor
end

@testset "non-homogeneous dimension solves the Moran equation" begin
    Γ = Attractor(Similarity(1 / 3, 0.0), Similarity(1 / 2, 1 / 2))
    @test sum(s.ρ^Γ.d for s in Γ.ifs) ≈ 1
    μ = HausdorffMeasure(Γ)
    @test μ.weights ≈ [s.ρ^Γ.d for s in Γ.ifs]
end

@testset "preset $name" for name in ALL_PRESETS
    Γ = getfractal(name)
    μ = HausdorffMeasure(Γ)
    M = length(Γ.ifs)

    @testset "Moran equation, and weights" begin
        @test sum(s.ρ^Γ.d for s in Γ.ifs) ≈ 1
        @test sum(μ.weights) ≈ 1
        @test all(μ.weights .≈ [s.ρ^Γ.d for s in Γ.ifs])
    end

    @testset "contractions" begin
        @test all(0 < s.ρ < 1 for s in Γ.ifs)
    end

    @testset "barycentre is the weighted barycentre of its children" begin
        @test F.get_barycentre(μ) ≈ sum(μ.weights[m] * F.get_barycentre(μ[m]) for m = 1:M) atol =
            1e-12
    end

    @testset "sub-measures" begin
        @test sum(μ[m].suppmeasure for m = 1:M) ≈ μ.suppmeasure
        @test μ[1, 1].suppmeasure ≈ μ.weights[1]^2 * μ.suppmeasure
        @test Γ[1].diam ≈ Γ.ifs[1].ρ * Γ.diam
        @test Γ[1, 2].diam ≈ Γ.ifs[1].ρ * Γ.ifs[2].ρ * Γ.diam
    end

    @testset "first mesh level partitions the measure" begin
        V = F.construct_quasiuniform_prebary_p0basis(μ, 0.4 * Γ.diam, 0.2 * Γ.diam)
        @test sum(ϕ.measure.suppmeasure for ϕ in V) ≈ μ.suppmeasure
        @test all(ϕ.measure.supp.diam < 0.4 * Γ.diam for ϕ in V)
        # L² normalised
        @test all(ϕ ⋅ ϕ ≈ 1 for ϕ in V)
    end
end

@testset "connectedness matrix of $name is symmetric" for name in (
    "cantorset",
    "cantordust",
    "sierpinskitriangle",
    "sierpinskicarpet",
    "vicsek",
    "kochcurve",
    "kochsnowflake",
    "heighwaydragon",
)
    C = getfractal(name).connectedness
    @test C == C'
end

@testset "symmetries of $name map the attractor to itself" for (name, kw) in (
    ("cantorset", (;)),
    ("cantordust", (;)),
    ("cantordust", (ρ = 1 / 4,)),
    ("sierpinskitriangle", (;)),
    ("sierpinskitriangle", (ρ = 1 / 3,)),
    ("sierpinskicarpet", (;)),
    ("vicsek", (;)),
    ("kochcurve", (;)),
    ("kochsnowflake", (;)),
)
    Γ = getfractal(name; kw...)
    h = name in ("kochsnowflake", "sierpinskicarpet") ? 0.3 : 0.1
    pts, _ = F.barycentre_quadrule(HausdorffMeasure(Γ), h * Γ.diam)
    # the rule's nodes are barycentres of a mesh, which symmetries permute
    @test maximum(minimum(norm(g(x) - y) for y in pts) for g in Γ.symmetries, x in pts) <
          1e-8
end

@testset "translation" begin
    Γ = getfractal("cantordust")
    v = [2.0, -1.0]
    Γᵥ = Γ + v
    @test Γᵥ.diam == Γ.diam
    @test Γᵥ.d == Γ.d
    @test F.get_barycentre(HausdorffMeasure(Γᵥ)) ≈ F.get_barycentre(HausdorffMeasure(Γ)) + v
    @test (v+Γ).ifs[2].δ ≈ Γᵥ.ifs[2].δ
    @test_throws AssertionError Γ + [1.0, 2.0, 3.0]

    Γ₁ = getfractal("cantorset") + 1.5
    @test F.get_barycentre(HausdorffMeasure(Γ₁)) ≈ 2
end

@testset "embedding" begin
    Γ = getfractal("cantorset")
    E = F.embed(Γ, 2)
    @test F.get_ambient_dimension(E) == 2
    @test E.d == Γ.d
    @test E.diam == Γ.diam
    @test F.embed(Γ, 1) === Γ
    x = [0.3, 0.0]
    @test E.ifs[2](x) ≈ [Γ.ifs[2](0.3), 0.0]
    @test_throws AssertionError F.embed(getfractal("cantordust"), 1)
    μ = HausdorffMeasure(Γ)
    @test F.get_barycentre(F.embed(μ, 3)) ≈ [0.5, 0.0, 0.0]
end

@testset "unions" begin
    Γ₁ = getfractal("cantordust")
    Γ₂ = getfractal("cantorset")
    U = Γ₁ ∪ Γ₂
    @test length(U) == 2
    @test F.get_ambient_dimension(U[2]) == 2
    @test length(U ∪ Γ₁) == 3

    μ₁ = HausdorffMeasure(Γ₁)
    μ₂ = HausdorffMeasure(Γ₂)
    V = μ₁ ∪ μ₂
    @test length(V) == 2
    @test F.get_ambient_dimension(V[2]) == 2
    # nothing to embed, so nothing is rebuilt
    @test (μ₁∪μ₁)[1] === μ₁

    @testset "bases on unions" begin
        V₁ = F.construct_quasiuniform_prebary_p0basis(μ₁, 0.4, 0.2)
        V₂ = F.construct_quasiuniform_prebary_p0basis(μ₁, 0.2, 0.1)
        B = V₁ ∪ V₂
        @test length(B) == length(V₁) + length(V₂)
        @test size(B) == (length(B),)
        @test B[1] === V₁[1]
        @test B[length(V₁)+1] === V₂[1]
        @test B[length(B)] === V₂[end]
        @test_throws AssertionError B[0]
        @test F.get_ambient_dimension(B) == 2
    end
end

@testset "vector indices" begin
    Γ = getfractal("cantordust")
    μ = HausdorffMeasure(Γ)

    @testset "round trip" for 𝐦 in ([1], [3, 1], [2, 4, 1, 3], [4, 4, 4, 4, 4])
        v = F.VectorIndex{4}(𝐦)
        @test Vector(v) == 𝐦
        @test length(v) == length(𝐦)
        @test [v[j] for j in eachindex(𝐦)] == 𝐦
    end

    @testset "mesh of $name" for name in ("cantorset", "cantordust", "kochcurve", "vicsek")
        Γₙ = getfractal(name)
        h = Γₙ.diam / 10
        L = F.subdivide_indices(Γₙ, h)
        @test all(Γₙ[Vector(𝐦)].diam < h for 𝐦 in L)
        # no element is inside another, and together they exhaust the measure
        @test length(unique(Vector.(L))) == length(L)
        μₙ = HausdorffMeasure(Γₙ)
        @test sum(μₙ[Vector(𝐦)].suppmeasure for 𝐦 in L) ≈ μₙ.suppmeasure
    end

    @testset "non-homogeneous mesh is still quasi-uniform" begin
        Γₙ = Attractor(Similarity(1 / 3, 0.0), Similarity(1 / 2, 1 / 2))
        h = 0.05
        L = F.subdivide_indices(Γₙ, h)
        @test all(Γₙ[Vector(𝐦)].diam < h for 𝐦 in L)
        @test minimum(Γₙ[Vector(𝐦)].diam for 𝐦 in L) > h / 3
    end
end
