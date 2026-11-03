# Quadrature on measures whose moments are known

@testset "barycentre rule" begin
    @testset "exact for linear functions on $name" for name in (
        "cantorset",
        "cantordust",
        "kochcurve",
    )
        Γ = getfractal(name)
        μ = HausdorffMeasure(Γ)
        x, w = F.barycentre_quadrule(μ, 0.05)
        @test sum(w) ≈ μ.suppmeasure
        @test sum(w .* x) ≈ F.get_barycentre(μ) * μ.suppmeasure
    end

    @testset "converges on the interval" begin
        Γ = getfractal("cantorset", ρ = 1 / 2)
        err(h) = begin
            x, w = F.barycentre_quadrule(Γ, h)
            return abs(sum(w .* x .^ 2) - 1 / 3)
        end
        @test err(0.01) < err(0.1) / 10
        @test err(0.01) < 1e-3
    end
end

@testset "Gauss rule" begin
    @testset "matches Gauss-Legendre on the interval" begin
        Γ = getfractal("cantorset", ρ = 1 / 2)
        x, w = F.gauss_quadrule(Γ, 3)
        @test sort(x) ≈ [(1 - sqrt(3 / 5)) / 2, 1 / 2, (1 + sqrt(3 / 5)) / 2]
        @test w[sortperm(x)] ≈ [5, 8, 5] / 18
    end

    @testset "exact for polynomials on the Cantor set" begin
        # moments from self-similarity: mₙ = (ρⁿ/2)Σ C(n,k) (1-ρ)^(n-k)/ρ^(n-k) mₖ + (ρⁿ/2) mₙ
        ρ = 1 / 3
        m = [1.0]
        for n = 1:9
            rhs = sum(binomial(n, k) * ρ^k * (1 - ρ)^(n - k) * m[k+1] for k = 0:(n-1))
            push!(m, rhs / (2 - 2 * ρ^n))
        end
        x, w = F.gauss_quadrule(getfractal("cantorset"), 5)
        @testset "degree $n" for n = 0:9
            @test sum(w .* x .^ n) ≈ m[n+1] rtol = 1e-10
        end
    end

    @test all(F.gauss_quadrule(getfractal("cantorset"), 6)[2] .> 0)
    @test_throws AssertionError F.gauss_quadrule(getfractal("cantordust"), 3)

    # A Jacobi matrix of the normalised measure is what's needed, then scale the weights.
    @testset "measures of mass other than one" begin
        μ = HausdorffMeasure(getfractal("cantorset"))
        @test_broken begin
            x, w = F.gauss_quadrule(μ[1], 3)
            sum(w) ≈ μ[1].suppmeasure
        end
    end
end

@testset "chaos game" begin
    Random.seed!(1234)
    μ = HausdorffMeasure(getfractal("cantordust"))
    x, w = F.chaos_quadrule(μ, 100_000)
    @test sum(w) ≈ μ.suppmeasure
    @test sum(w .* x) ≈ F.get_barycentre(μ) atol = 0.01
end

@testset "tensor quadrature" begin
    μ = HausdorffMeasure(getfractal("cantorset"))
    x₁, w₁ = F.gauss_quadrule(μ, 3)
    x, y, w = F.combine_quadrules(x₁, w₁, x₁, w₁)
    @test length(x) == length(y) == length(w) == 9
    @test sum(w) ≈ 1
    @test sum(w .* x .* y) ≈ 1 / 4
end

@testset "s-energy" begin
    @testset "exact on the interval" begin
        # ∫∫log|x-y| = -3/2, and ∫∫|x-y|^-s = 2/((1-s)(2-s))
        Γ = getfractal("cantorset", ρ = 1 / 2)
        @test F.s_energy(Γ, 0.0; h_quad = 0.01) ≈ -3 / 2 rtol = 1e-3
        for s in (0.25, 0.5, 0.75)
            @test F.s_energy(Γ, s; h_quad = 0.01) ≈ 2 / ((1 - s) * (2 - s)) rtol = 1e-3
        end
    end

    @testset "exact on the square" begin
        # mean inverse distance between points of the unit square
        Γ = getfractal("cantordust", ρ = 1 / 2)
        exact = 4 * (log(1 + sqrt(2)) + (1 - sqrt(2)) / 3)
        @test F.s_energy(Γ, 1.0; h_quad = 0.05) ≈ exact rtol = 1e-2
    end

    @testset "scaling of the energy on $name" for name in ("cantorset", "cantordust")
        # rescaling Γ by λ multiplies the energy by λ^-s
        Γ = getfractal(name)
        μ = HausdorffMeasure(Γ)
        E = F.s_energy(μ, 0.5; h_quad = 0.02)
        @test F.s_energy(μ[1], 0.5; h_quad = 0.02 * Γ.ifs[1].ρ) ≈
              E * μ.weights[1]^2 * Γ.ifs[1].ρ^(-0.5) rtol = 1e-8
    end

    @testset "self-similar log relation on the Cantor set" begin
        # the log kernel gains a constant under rescaling, rather than a factor
        Γ = getfractal("cantorset")
        μ = HausdorffMeasure(Γ)
        E = F.s_energy(μ, 0.0; h_quad = 0.02)
        E₁ = F.s_energy(μ[1], 0.0; h_quad = 0.02 * Γ.ifs[1].ρ)
        @test E₁ ≈ μ.weights[1]^2 * (E + log(Γ.ifs[1].ρ)) rtol = 1e-8
    end

    @testset "symmetry only reduces the work, on $name" for name in
                                                            ("cantorset", "kochcurve")
        Γ = getfractal(name)
        Γₜ = Attractor(
            collect(Γ.ifs);
            diam = Γ.diam,
            d = Γ.d,
            connectedness = Γ.connectedness,
            symmetries = F.trivialgroup(eltype(Γ.ifs[1].ρ), F.get_ambient_dimension(Γ)),
        )
        @test F.s_energy(Γ, 0.5; h_quad = 0.02) ≈ F.s_energy(Γₜ, 0.5; h_quad = 0.02) rtol =
            1e-8
    end
end
