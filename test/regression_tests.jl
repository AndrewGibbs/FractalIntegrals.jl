# One test set per fixed bug: see dev_notes/changes_to_existing_code.md
import LinearAlgebra: dot, norm
const F = FractalIntegrals

@testset "block operators: ∪ inside the package" begin
    Γ₁ = getfractal("cantor dust")
    Γ₂ = getfractal("cantor set") + 1.5
    S = F.singlelayer_operator_helmholtz((Γ₁, Γ₂), 2.0)
    @test S isa F.BlockOperator
    @test F.singlelayer_operator_laplace((Γ₁, Γ₂)) isa F.BlockOperator
    bases = (
        F.construct_quasiuniform_prebary_p0basis(S.measures[1], 0.5, 0.1),
        F.construct_quasiuniform_prebary_p0basis(S.measures[2], 0.5, 0.1),
    )
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
    S₀ = F.singlelayer_operator_helmholtz(Γ, 1e-6; ambient_dimension = 3)
    L = F.singlelayer_operator_laplace(Γ; ambient_dimension = 3)
    A₀ = F.discretise(S₀; h_mesh = h, h_quad = h_quad).stiffness_matrix
    Aₗ = F.discretise(L; h_mesh = h, h_quad = h_quad).stiffness_matrix
    @test norm(A₀ - Aₗ) / norm(Aₗ) < 1e-5

    k = 3.0
    S = F.singlelayer_operator_helmholtz(Γ, k; ambient_dimension = 3)
    Sₕ = F.discretise(S; h_mesh = h, h_quad = h_quad)

    # near field must tend to far field as r → ∞
    ϕ = Sₕ \ (x -> exp(im * k * x[1]))
    slp = F.singlelayer_potential_helmholtz(ϕ, k; ambient_dimension = 3)
    ffp = F.farfield_pattern(ϕ, k; ambient_dimension = 3)
    r = 500
    angles = [(θ, ψ) for θ = 0.3:0.5:2.8, ψ = 0.0:1.0:5.0][:]
    x = [r * [sin(θ) * cos(ψ), sin(θ) * sin(ψ), cos(θ)] for (θ, ψ) in angles]
    @test (exp(im * k * r) / r) * ffp.(angles) ≈ -slp.(x) rtol = 1e-2

    # reciprocity of the far-field pattern
    dirs = [(0.4, 0.3), (0.9, 1.7), (1.3, 3.1), (0.7, 4.6)]
    fθ = [x -> exp(-im * k * sin(θ) * (x[1] * cos(ψ) + x[2] * sin(ψ))) for (θ, ψ) in dirs]
    R = [
        F.farfield_pattern(Sₕ \ fθ[n], k; ambient_dimension = 3)(dirs[m]) for
        m in eachindex(dirs), n in eachindex(dirs)
    ]
    @test R ≈ transpose(R) rtol = 1e-2

    # also for a fractal with ambient dimension three
    cube = Attractor(
        [
            Similarity(1 / 3, [a, b, c]) for a in (0, 2 / 3), b in (0, 2 / 3),
            c in (0, 2 / 3)
        ][:]...,
    )
    @test F.singlelayer_operator_helmholtz(cube, k) isa
          F.OscillatorySeparableIntegralOperator
end

@testset "Hausdorff dimension of presets" begin
    Γ = getfractal("sierpinski")
    @test Γ.d ≈ log(3) / log(2)
    @test sum(HausdorffMeasure(Γ).weights) ≈ 1
    @test length(F.chaos_quadrule(Γ, 1000)[1]) == 1000
    @test F.sierpinskicarpet(; ρ = 1 / 4).d ≈ 3 / 2
    # every preset without a user-set dimension rule must have weights summing to one
    for name in ["cantorset", "cantordust", "sierpinski", "vicsek", "kochcurve", "carpet"]
        @test sum(HausdorffMeasure(getfractal(name)).weights) ≈ 1
    end
end

@testset "barycentre rule, non-Hausdorff measure on homogeneous attractor" begin
    μ = InvariantMeasure(getfractal("cantordust"), [0.1, 0.2, 0.3, 0.4])
    x, w = F.barycentre_quadrule(μ, 0.1)
    @test length(x) == length(w)
    @test sum(w) ≈ 1
    # the first moment is integrated exactly, which also checks weights match their nodes
    @test sum(w .* x) ≈ F.get_barycentre(μ)
end

@testset "getmeasure with unknown measure" begin
    μ = @test_logs (:warn, r"not known") getmeasure("cantordust")
    @test μ isa HausdorffMeasure
    @test μ.suppmeasure == 1
    @test getmeasure("cantordust"; suppmeasure = 2.0).suppmeasure == 2
end

@testset "single-element mesh" begin
    Γ = getfractal("cantordust")
    @test F.discretise(F.IdentityOperator(Γ); h_mesh = 10.0).stiffness_matrix ≈ ones(1, 1)
end

@testset "inner product of elements at different depths" begin
    μ = HausdorffMeasure(getfractal("cantordust"))
    coarse = F.construct_quasiuniform_prebary_p0basis(μ, 0.3, 0.1)
    fine = F.construct_quasiuniform_prebary_p0basis(μ, 0.1, 0.05)
    a = coarse[1]
    inside =
        first(e for e in fine if Vector(e.vindex)[1:length(a.vindex)] == Vector(a.vindex))
    outside =
        first(e for e in fine if Vector(e.vindex)[1:length(a.vindex)] != Vector(a.vindex))
    # the overlap is the smaller element
    expected = inside.measure.suppmeasure * a.normalisation * inside.normalisation
    @test dot(a, inside) ≈ expected
    @test dot(inside, a) ≈ expected
    @test dot(a, outside) == 0
end

@testset "small definition errors" begin
    Γ = getfractal("cantordust")
    μ = HausdorffMeasure(Γ)
    @test length(Base.:∪(μ ∪ μ, μ)) == 3
    s = Similarity(1 / 3, [0.2, 0.1], [0.0 -1.0; 1.0 0.0])
    @test F.sim_map_inv(s, s([0.3, 0.7])) ≈ [0.3, 0.7]
    @test all(n -> isdefined(F, n), names(F))
end

@testset "getindex only extended for BlockOperator" begin
    @test_throws MethodError getindex(:not_indexable, 1, 2)
    Γ = getfractal("cantordust")
    S = F.singlelayer_operator_laplace((Γ, Γ + [3.0, 0.0]))
    @test S[1, 2] === S.operators[1, 2]
end

@testset "docstring example for Attractor" begin
    courage = Similarity(1 / 2, [0, 0])
    wisdom = Similarity(1 / 2, [1 / 2, 0])
    power = Similarity(1 / 2, [1 / 4, sqrt(3) / 4])
    Γ = Attractor(courage, wisdom, power; connectedness = ones(Bool, 3, 3))
    @test Γ.d ≈ log(3) / log(2)
end

@testset "block operators solve and give fields" begin
    f(x) = 1.0 + 0im
    # Splitting a fractal into its components should change nothing but the block layout
    relerr(A, B) = norm(A - B) / norm(A)
    matrix(K; kw...) = F.discretise(K; h_mesh = 0.1, h_quad = 0.03, kw...).stiffness_matrix

    @testset "Cantor set in two parts" begin
        μ = HausdorffMeasure(getfractal("cantorset"))
        pieces = (μ[1], μ[2])
        @test relerr(
            matrix(F.singlelayer_operator_laplace(μ)),
            matrix(F.singlelayer_operator_laplace(pieces)),
        ) < 1e-13
        @test relerr(
            matrix(F.singlelayer_operator_helmholtz(μ, 5.0)),
            matrix(F.singlelayer_operator_helmholtz(pieces, 5.0)),
        ) < 1e-13
    end

    @testset "Cantor dust in four parts" begin
        μ = HausdorffMeasure(getfractal("cantordust"))
        pieces = (μ[1], μ[2], μ[3], μ[4])
        K = F.singlelayer_operator_laplace(pieces)
        @test relerr(matrix(F.singlelayer_operator_laplace(μ)), matrix(K)) < 1e-13
        @test eltype(matrix(K)) <: Real
        @test matrix(K) ≈ matrix(K; reps = false)

        # the choice of ambient dimension must reach the diagonal blocks too
        K₃ = F.singlelayer_operator_laplace(pieces, ambient_dimension = 3)
        @test all(K₃[n, n].s == 1 for n = 1:4)
        S₃ = F.singlelayer_operator_laplace(μ, ambient_dimension = 3)
        @test relerr(matrix(S₃), matrix(K₃)) < 1e-13
    end

    @testset "solution and fields" begin
        μ = HausdorffMeasure(getfractal("cantordust"))
        pieces = (μ[1], μ[2], μ[3], μ[4])
        k = 5.0
        ϕ = F.discretise(F.singlelayer_operator_helmholtz(μ, k), h_mesh = 0.1, h_quad = 0.03) \ f
        ϕ₄ = F.discretise(
            F.singlelayer_operator_helmholtz(pieces, k),
            h_mesh = 0.1,
            h_quad = 0.03,
        ) \ f
        x = [0.3, 1.2]
        @test F.singlelayer_potential_helmholtz(ϕ₄, k)(x) ≈
              F.singlelayer_potential_helmholtz(ϕ, k)(x)
        @test F.farfield_pattern(ϕ₄, k)(0.3) ≈ F.farfield_pattern(ϕ, k)(0.3)
    end

    @testset "default mesh width, and one-dimensional sets with ambient dimension one" begin
        Γ = getfractal("cantorset")
        K = F.singlelayer_operator_laplace((Γ, Γ + 3.0), ambient_dimension = 1)
        ϕ = F.discretise(K) \ f
        @test isfinite(F.singlelayer_potential_laplace(ϕ)([0.5, 1.0]))
    end

    @testset "union of sets with different ambient dimensions" begin
        # The separation is chosen so that the log capacity of the union is far from one.
        # (At a separation of 2 the block matrix is nearly singular, as it should be.)
        Γ₁ = getfractal("cantorset")
        Γ₂ = getfractal("cantordust") + [0.0, 1.2]
        K = F.singlelayer_operator_laplace((Γ₁, Γ₂))
        @test all(F.get_ambient_dimension(K.measures[n]) == 2 for n = 1:2)
        ϕ = F.discretise(K, h_mesh = 0.05, h_quad = 0.0125) \ f
        slp = F.singlelayer_potential_laplace(ϕ)
        @test real(slp([0.0, 1.2])) ≈ 1 atol = 0.1 # corner of the dust
        @test real(slp([1 / 3, 0.0])) ≈ 1 atol = 0.2 # a point of the Cantor set

        # embedding must not change a matrix
        S₁ = F.singlelayer_operator_laplace(HausdorffMeasure(Γ₁))
        S₂ = F.singlelayer_operator_laplace(HausdorffMeasure(F.embed(Γ₁, 2)))
        @test matrix(S₁) ≈ matrix(S₂)
    end

    @testset "collocation is refused with a clear error" begin
        μ = HausdorffMeasure(getfractal("cantorset"))
        K = F.singlelayer_operator_laplace((μ[1], μ[2]))
        @test_throws ArgumentError F.discretise(K, method = :collocation)
    end
end

@testset "SingularIntegralOperator" begin
    Γ = getfractal("cantorset")
    μ = HausdorffMeasure(Γ)
    Φ(x, y) = -F.energykernel(0, x, y) / (2π)
    K = F.SingularIntegralOperator(μ, Φ)
    @test K isa F.IntegralOperator
    @test F.SingularIntegralOperator(Γ, Φ, true).symmetric

    # must match the same kernel, written as the Laplace operator
    S = F.singlelayer_operator_laplace(μ)
    Kₕ = F.discretise(K, h_mesh = 0.1, h_quad = 0.01, method = :collocation)
    Sₕ = F.discretise(S, h_mesh = 0.1, h_quad = 0.01, method = :collocation)
    @test Kₕ.collocation_matrix ≈ Sₕ.collocation_matrix

    @test_throws ArgumentError F.discretise(K, h_mesh = 0.1, h_quad = 0.01)
end

@testset "second-kind equations under collocation" begin
    Γ = getfractal("cantordust")
    μ = HausdorffMeasure(Γ)
    I = F.IdentityOperator(μ)
    S = F.singlelayer_operator_laplace(μ)
    opts = (h_mesh = 0.1, h_quad = 0.01, method = :collocation)

    @testset "identity" begin
        Iₕ = F.discretise(I; opts...)
        Vₕ = Iₕ.basis
        @test Iₕ.collocation_matrix ≈ F.Diagonal([ϕ.normalisation for ϕ in Vₕ])
        # the whole attractor, index [0], contains every collocation point
        Mₕ = F.get_collocation_matrix(
            I,
            F.construct_quasiuniform_prebary_p0basis(μ, 10.0, 0.5),
            0.1,
            0.01,
        )[1]
        @test all(Mₕ .== 1)
    end

    @testset "linearity" begin
        Sₕ = F.discretise(S; opts...).collocation_matrix
        Iₕ = F.discretise(I; opts...).collocation_matrix
        Aₕ = F.discretise(2.5 * I + S; opts...).collocation_matrix
        @test Aₕ ≈ 2.5 * Iₕ + Sₕ
    end

    @testset "agrees with Galerkin" begin
        # entries differ, but they solve the same equation
        f(x) = 1.0 + x[1] + 0im
        A = 2 * I + S
        c = F.discretise(A; opts...) \ f
        g = F.discretise(A, h_mesh = 0.1, h_quad = 0.01) \ f
        @test c.coeffs ≈ g.coeffs rtol = 0.05
    end

    @testset "keywords and default mesh widths" begin
        A = 2 * I + S
        @test F.getdefault_meshwidth(A) == F.getdefault_meshwidth(S)
        Mᵣ = F.discretise(A, h_mesh = 0.1, h_quad = 0.01, reps = true).stiffness_matrix
        M = F.discretise(A, h_mesh = 0.1, h_quad = 0.01, reps = false).stiffness_matrix
        @test Mᵣ ≈ M
    end
end
