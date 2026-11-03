# One test set per fixed bug: see dev_notes/changes_to_existing_code.md
import LinearAlgebra: dot
const F = FractalIntegrals

@testset "block operators: ∪ inside the package" begin
    Γ₁ = getfractal("cantor dust")
    Γ₂ = getfractal("cantor set") + 1.5
    S = F.singlelayer_operator_helmholtz((Γ₁, Γ₂), 2.0)
    @test S isa F.BlockOperator
    @test F.singlelayer_operator_laplace((Γ₁, Γ₂)) isa F.BlockOperator
    bases = (F.construct_quasiuniform_prebary_p0basis(S.measures[1], 0.5, 0.1),
            F.construct_quasiuniform_prebary_p0basis(S.measures[2], 0.5, 0.1))
    G = F.get_galerkin_matrix(S, bases)
    n₁ = length(bases[1])
    @test size(G) == (sum(length.(bases)), sum(length.(bases)))
    @test G[1:n₁, 1:n₁] ≈ F.get_galerkin_matrix(S[1, 1], bases[1])
    @test G ≈ transpose(G)
    # the polygon method of ∪ must still work
    @test length(F.polygon_approx(getfractal("kochsnowflake")).nodes) > 0
end

@testset "3D Helmholtz operator" begin
    Γ = getfractal("cantordust")
    h = 0.3
    h_quad = 0.05

    # as k → 0 the Helmholtz matrix must tend to the Laplace matrix
    S₀ = F.singlelayer_operator_helmholtz(Γ, 1e-6, ambient_dimension = 3)
    L = F.singlelayer_operator_laplace(Γ, ambient_dimension = 3)
    A₀ = F.discretise(S₀, h_mesh = h, h_quad = h_quad).stiffness_matrix
    Aₗ = F.discretise(L, h_mesh = h, h_quad = h_quad).stiffness_matrix
    @test norm(A₀ - Aₗ) / norm(Aₗ) < 1e-5

    k = 3.0
    S = F.singlelayer_operator_helmholtz(Γ, k, ambient_dimension = 3)
    Sₕ = F.discretise(S, h_mesh = h, h_quad = h_quad)

    # near field must tend to far field as r → ∞
    ϕ = Sₕ \ (x -> exp(im * k * x[1]))
    slp = F.singlelayer_potential_helmholtz(ϕ, k, ambient_dimension = 3)
    ffp = F.farfield_pattern(ϕ, k, ambient_dimension = 3)
    r = 500
    angles = [(θ, ψ) for θ in 0.3:0.5:2.8, ψ in 0.0:1.0:5.0][:]
    x = [r * [sin(θ) * cos(ψ), sin(θ) * sin(ψ), cos(θ)] for (θ, ψ) in angles]
    @test (exp(im * k * r) / r) * ffp.(angles) ≈ -slp.(x) rtol = 1e-2

    # reciprocity of the far-field pattern
    dirs = [(0.4, 0.3), (0.9, 1.7), (1.3, 3.1), (0.7, 4.6)]
    fθ = [x -> exp(-im * k * sin(θ) * (x[1] * cos(ψ) + x[2] * sin(ψ))) for (θ, ψ) in dirs]
    R = [F.farfield_pattern(Sₕ \ fθ[n], k, ambient_dimension = 3)(dirs[m]) for m in eachindex(dirs), n in eachindex(dirs)]
    @test R ≈ transpose(R) rtol = 1e-2

    # also for a fractal with ambient dimension three
    cube = Attractor([Similarity(1/3, [a, b, c]) for a in (0, 2/3), b in (0, 2/3), c in (0, 2/3)][:]...)
    @test F.singlelayer_operator_helmholtz(cube, k) isa F.OscillatorySeparableIntegralOperator
end

@testset "Hausdorff dimension of presets" begin
    Γ = getfractal("sierpinski")
    @test Γ.d ≈ log(3) / log(2)
    @test sum(HausdorffMeasure(Γ).weights) ≈ 1
    @test length(F.chaos_quadrule(Γ, 1000)[1]) == 1000
    @test F.sierpinskicarpet(ρ = 1/4).d ≈ 3/2
    # every preset without a user-set dimension rule must have weights summing to one
    for name in ["cantorset", "cantordust", "sierpinski", "vicsek", "kochcurve", "carpet"]
        @test sum(HausdorffMeasure(getfractal(name)).weights) ≈ 1
    end
end
