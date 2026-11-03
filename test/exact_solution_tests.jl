# Checks against answers known independently of this package.
# Each test treats a non-fractal shape (a segment, or a square) as a self-similar set.

@testset "Laplace on a segment in ℝ², exact solution" begin
    # Solving Sφ = 1 on [0,1]×{0}: the exact field comes from the Joukowski map,
    # u(z) = m⋅U(z), where m = 2π/log(4) is the total mass of the density φ.
    Γ = Attractor(
        Similarity(1 / 2, [0.0, 0.0]),
        Similarity(1 / 2, [0.5, 0.0]);
        connectedness = ones(Bool, 2, 2),
    )
    function U(p)
        ζ = 2 * complex(p[1], p[2]) - 1
        return -(log(abs((ζ + sqrt(ζ - 1) * sqrt(ζ + 1)) / 2)) - log(2)) / (2π)
    end
    m = 2π / log(4)
    S = F.singlelayer_operator_laplace(Γ)
    ϕ = F.discretise(S, h_mesh = 0.01, h_quad = 0.0025) \ (x -> 1.0 + 0im)

    @testset "total mass of density" begin
        @test real(density_mass(ϕ)) ≈ m rtol = 1e-2
    end

    slp = F.singlelayer_potential_laplace(ϕ)
    points = ([0.5, 1.0], [1.5, 0.5], [-1.0, 2.0], [3.0, 3.0], [0.5, 0.05])
    @testset "field at $p" for p in points
        @test real(slp(p)) ≈ m * U(p) atol = 2e-3
    end

    @testset "field tends to boundary data on the segment" begin
        @test real(slp([0.5, 1e-3])) ≈ 1 atol = 2e-2
        @test abs(imag(slp([0.5, 1.0]))) < 1e-12
    end

    # The same problem, posed on a fractal with ambient dimension one
    Γ₁ = getfractal("cantorset", ρ = 1 / 2)
    ϕ₁ = F.discretise(F.singlelayer_operator_laplace(Γ₁), h_mesh = 0.01, h_quad = 0.0025) \
         (x -> 1.0 + 0im)
    slp₁ = F.singlelayer_potential_laplace(ϕ₁)
    @testset "screen potential, ambient dimension one, at $p" for p in points
        @test real(slp₁(p)) ≈ m * U(p) atol = 2e-3
    end
end

@testset "Laplace on the unit square as a screen in ℝ³" begin
    # The capacitance of a unit square plate is known to be about 0.3667.
    Γ = getfractal("cantordust", ρ = 1 / 2)
    S = F.singlelayer_operator_laplace(Γ, ambient_dimension = 3)
    capacitance(h) = begin
        ϕ = F.discretise(S, h_mesh = h, h_quad = h / 4) \ (x -> 1.0 + 0im)
        real(density_mass(ϕ)) / (4π), ϕ
    end
    C₁, _ = capacitance(0.1)
    C₂, ϕ = capacitance(0.05)
    # error is observed to halve with h, so extrapolate
    @test 2C₂ - C₁ ≈ 0.3667 rtol = 2e-3

    # far from the plate, the field must look like C/r
    slp = F.singlelayer_potential_laplace(ϕ, ambient_dimension = 3)
    r = 2000.0
    x = [0.5, 0.5, 0.0] + r * [0.3, -0.5, 0.8] / norm([0.3, -0.5, 0.8])
    @test r * real(slp(x)) ≈ C₂ rtol = 1e-3
    # and near the plate, it must tend to the boundary data
    @test real(slp([0.5, 0.5, 1e-3])) ≈ 1 atol = 2e-2
end

@testset "Laplace potential of a density on the unit square in ℝ²" begin
    # Since u = Sφ is a Newton potential, -Δu = φ inside the square and Δu = 0 outside
    Γ = getfractal("cantordust", ρ = 1 / 2)
    μ = HausdorffMeasure(Γ)
    V = F.construct_quasiuniform_prebary_p0basis(μ, 0.03, 0.015)
    φ(x) = sin(π * x[1]) * sin(π * x[2])
    u = F.singlelayer_potential_laplace(F.project(V, x -> φ(x) + 0im))
    H = 0.15
    # five-point stencil, whose error for φ is a factor of about 1-π²H²/12 here
    Δₕ(p) = real(u(p + [H, 0]) + u(p - [H, 0]) + u(p + [0, H]) + u(p - [0, H]) - 4u(p)) / H^2

    @testset "-Δu = φ at interior point $p" for p in ([0.37, 0.41], [0.5, 0.5], [0.7, 0.25])
        @test -Δₕ(p) ≈ φ(p) rtol = 5e-2
    end
    @testset "Δu = 0 at exterior point $p" for p in ([1.6, 0.5], [0.5, -0.8])
        @test abs(Δₕ(p)) < 5e-3
    end

    # logarithmic growth, with the scaling of the free-space Green's function
    u₁ = F.singlelayer_potential_laplace(F.project(V, x -> 1.0 + 0im))
    c = [0.5, 0.5]
    @test real(u₁(c + [10.0, 0]) - u₁(c + [100.0, 0])) ≈ log(10) / (2π) rtol = 1e-3

    # the k → 0 limit of the Helmholtz potential differs only by a constant
    hel = F.singlelayer_potential_helmholtz(F.project(V, x -> 1.0 + 0im), 1e-6)
    p, q = [2.0, 0.0], [0.5, 3.0]
    @test real(hel(p) - hel(q)) ≈ real(u₁(p) - u₁(q)) rtol = 1e-3
end
