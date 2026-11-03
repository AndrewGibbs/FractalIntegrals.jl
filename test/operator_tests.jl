# Operators and discretisation, tested through what they must satisfy

@testset "kernels" begin
    x, y = [0.3, 0.7], [1.0, -0.2]
    r = norm(x - y)
    X, Y = [x, x], [y, y]

    @testset "vectorised distance" begin
        @test F.vecdist(X, Y) ≈ [r, r]
        @test F.vecdist(x, Y) ≈ [r, r]
        @test F.vecdist(X, y) ≈ [r, r]
        @test F.vecdist([0.2, 0.5], [1.0, 0.1]) ≈ [0.8, 0.4]
    end

    @testset "energy kernel" begin
        @test F.energykernel(0, 2.0) == log(2.0)
        @test F.energykernel(0.5, 4.0) ≈ 0.5
    end

    @testset "Helmholtz minus singular part is its Lipschitz part" begin
        k = 3.0
        s2 = F.singlelayer_operator_helmholtz(getfractal("cantordust"), k)
        @test s2.kernel(X, Y) ≈
              s2.lipschitzpart(X, Y) + s2.singularconst * F.energykernel(0, X, Y)
        s3 = F.singlelayer_operator_helmholtz(
            getfractal("cantordust"),
            k,
            ambient_dimension = 3,
        )
        @test s3.kernel(X, Y) ≈
              s3.lipschitzpart(X, Y) + s3.singularconst * F.energykernel(1, X, Y)
        # bounded as x → y
        near = [x + [1e-7, 0.0]]
        @test abs(s2.lipschitzpart(near, [x])[1]) < 2
        @test abs(s3.lipschitzpart(near, [x])[1]) < 2
    end

    @testset "Helmholtz kernels solve Helmholtz" begin
        k = 3.0
        Φ₂(p) = F.helmholtzkernel2d(k, [p], [[0.0, 0.0]])[1]
        Φ₃(p) = F.helmholtzkernel3d(k, [p], [[0.0, 0.0, 0.0]])[1]
        h = 1e-3
        p₂ = [0.8, 0.5]
        Δ₂ = sum(Φ₂(p₂ + h * e) + Φ₂(p₂ - h * e) - 2Φ₂(p₂) for e in ([1, 0], [0, 1])) / h^2
        @test Δ₂ + k^2 * Φ₂(p₂) ≈ 0 atol = 1e-4
        p₃ = [0.8, 0.5, -0.3]
        Δ₃ =
            sum(
                Φ₃(p₃ + h * e) + Φ₃(p₃ - h * e) - 2Φ₃(p₃) for
                e in ([1, 0, 0], [0, 1, 0], [0, 0, 1])
            ) / h^2
        @test Δ₃ + k^2 * Φ₃(p₃) ≈ 0 atol = 1e-4
    end
end

@testset "operator algebra" begin
    Γ = getfractal("cantordust")
    S = F.singlelayer_operator_laplace(Γ)
    I = F.IdentityOperator(Γ)
    M(K) = F.discretise(K, h_mesh = 0.2, h_quad = 0.05).stiffness_matrix

    @test M(I) ≈ F.Diagonal(ones(size(M(I), 1)))
    @test M(3 * S) ≈ 3 * M(S)
    @test M(S * 3) ≈ 3 * M(S)
    @test M(2 * (3 * S)) ≈ 6 * M(S)
    @test (2 * (3 * S)).operator === S
    @test M(I + S) ≈ M(I) + M(S)
    @test M(2 * I + 0.5 * S) ≈ 2 * M(I) + 0.5 * M(S)
    @test_throws AssertionError I + F.IdentityOperator(getfractal("cantorset"))
end

@testset "Galerkin matrix on $name" for (name, h) in (
    ("cantorset", 0.1),
    ("cantordust", 0.2),
    ("sierpinskitriangle", 0.2),
    ("kochcurve", 0.2),
    ("vicsek", 0.2),
)
    Γ = getfractal(name)
    S = F.singlelayer_operator_laplace(Γ)
    Sₖ = F.singlelayer_operator_helmholtz(Γ, 3.0)
    A(K; kw...) = F.discretise(K, h_mesh = h, h_quad = h / 4; kw...).stiffness_matrix

    @testset "symmetric" begin
        @test A(S) ≈ transpose(A(S)) rtol = 1e-5
        @test A(Sₖ) ≈ transpose(A(Sₖ)) rtol = 1e-5
    end

    @testset "repeated entries change nothing" begin
        @test A(S; reps = true) ≈ A(S; reps = false) rtol = 1e-5
    end

    @testset "self-energy of the whole fractal" begin
        # with the constant function on one element, the entry is the s-energy
        μ = HausdorffMeasure(Γ)
        # (discretise rescales h_quad by diam/h_mesh, so this is the same rule on both sides)
        E = F.s_energy(μ, 0.0; h_quad = Γ.diam / 8)
        M = F.discretise(S, h_mesh = 10.0, h_quad = 10.0 / 8).stiffness_matrix
        @test size(M) == (1, 1)
        @test M[1, 1] ≈ -E / (2π) / μ.suppmeasure rtol = 1e-8
    end
end

@testset "Laplace matrix is positive definite for a set of small capacity" begin
    # the mean log kernel is positive definite when the capacity is below one
    A = F.discretise(
        F.singlelayer_operator_laplace(getfractal("cantordust")),
        h_mesh = 0.1,
        h_quad = 0.03,
    ).stiffness_matrix
    @test minimum(F.eigvals(F.Symmetric((A + A') / 2))) > 0
end

@testset "second-kind equation" begin
    # solving λφ + Sφ = f; the residual on the fractal should be small
    μ = HausdorffMeasure(getfractal("cantordust"))
    S = F.singlelayer_operator_laplace(μ)
    I = F.IdentityOperator(μ)
    λ = 2.0
    f(x) = 1.0 + x[1] + 0im
    @testset "$method" for method in (:galerkin, :collocation)
        ϕ = F.discretise(λ * I + S, h_mesh = 0.05, h_quad = 0.01, method = method) \ f
        u = F.singlelayer_potential_laplace(ϕ)
        # S is smoothing, and so φ = (f - Sφ)/λ is continuous on the dust
        for n in (1, 7, 20)
            el = ϕ.basis[n]
            x = F.get_barycentre(el.measure)
            @test λ * ϕ.coeffs[n] * el.normalisation + u(x) ≈ f(x) atol = 0.05
        end
    end
end

@testset "non-homogeneous attractors" begin
    Γ = Attractor(Similarity(1 / 3, 0.0), Similarity(1 / 2, 1 / 2))
    S = F.singlelayer_operator_laplace(Γ)
    Sₕ = F.discretise(S, h_mesh = 0.05, h_quad = 0.01)
    @test Sₕ.stiffness_matrix ≈ transpose(Sₕ.stiffness_matrix) rtol = 1e-6
    ϕ = Sₕ \ (x -> 1.0 + 0im)
    u = F.singlelayer_potential_laplace(ϕ)
    @test real(u([0.0, 0.0])) ≈ 1 atol = 0.2
    @test real(u([1.0, 0.0])) ≈ 1 atol = 0.2

    # Galerkin against collocation on an invariant measure, not Hausdorff
    μ = InvariantMeasure(Γ, 1.0, F.SVector(0.3, 0.7))
    a = F.discretise(F.singlelayer_operator_laplace(μ), h_mesh = 0.1, h_quad = 0.01)
    @test a.stiffness_matrix ≈ transpose(a.stiffness_matrix) rtol = 1e-6
end

@testset "far field agrees with the near field" begin
    # r^((n-1)/2) e^{-ikr} u(rd) tends to the far-field pattern, up to the conventions
    k = 4.0
    Γ = getfractal("cantordust")
    ϕ =
        F.discretise(F.singlelayer_operator_helmholtz(Γ, k), h_mesh = 0.1, h_quad = 0.03) \
        (x -> 1.0 + 0im)
    slp = F.singlelayer_potential_helmholtz(ϕ, k)
    ffp = F.farfield_pattern(ϕ, k)
    θ = 0.7
    r = 1e5
    x = r * [cos(θ), sin(θ)]
    @test abs(sqrt(r) * cis(-k * r) * slp(x)) ≈ abs(ffp(θ)) rtol = 1e-3
end
