using FractalIntegrals
using Test, MAT
import LinearAlgebra: norm, transpose

# test that fractal presets all load
@testset "Preset fractals" begin
    for Γ in keys(FractalIntegrals.fractaldict)
        strΓ = String(Γ)
        @testset "$strΓ" begin
            @test_nowarn FractalIntegrals.getfractal(Γ)
            @test_nowarn FractalIntegrals.getfractal(BigFloat, Γ)
        end
    end
end

# Approximate ∫₀¹∫₀¹ Φₛ(x,y) dx dy using s-energy on cantor set with ρ=1/2.
@testset "singular line segment" begin
    Γ = FractalIntegrals.cantorset(ρ = 1/2)
    for s in rand(5)
        I = 2/((1 - s)*(2 - s))
        @testset "s=$s" begin
            @test FractalIntegrals.s_energy(Γ, s, N_quad = 10) ≈ I
            @test FractalIntegrals.s_energy(Γ, s, h_quad = 0.001) ≈ I rtol=1e-4
        end
    end
end

# Andrea Moiola's 'prefractal BEM' data, produced using a different method
@testset "Comparison against prefractal BEM" begin
    # physical parameters
    k = 30.0
    d = [0.5000, -0.8660]
    Γ = getfractal("cantorset")
    
    # integral operator
    Sₖ = FractalIntegrals.singlelayer_operator_helmholtz(Γ, k, ambient_dimension = 2)

    # RHS
    f(x) = exp(im*k*d[1]*x)

    # points to test the ffp
    h_θ = 2π/300
    θ = 0:h_θ:2π
    pf_data = matread("Lebesgue_FF_vals.mat")
    pf_vals = pf_data["FF_vals"]

    # constant estimate
    C = 0.3

    # test for a range of meshwidths
    for h in 2.0 .^((-6:-1:-12))
        Sₖₕ = FractalIntegrals.discretise(Sₖ, h_mesh = h)
        ϕₕ = Sₖₕ\f
        ffp = FractalIntegrals.farfield_pattern(ϕₕ, k, ambient_dimension = 2)
        @testset "h=$h" begin
            @test norm(ffp.(θ) - pf_vals, Inf) / norm(pf_vals, Inf) < C*h^(Γ.d)
        end
    end
end

@testset "2D reciprocity test" begin
    for Γname in ["cantordust", "sierpinski"]
        for k in [1, 2, 5]
            Γ = getfractal(Γname)

            N = 10 # number of angles to test reciprocity
            recipangles = rand(N)

            # create RHS data for each incident angle
            fθ = [x -> exp(-im*k*(x[1]*cos(θ)+x[2]*sin(θ))) for θ in recipangles]

            # create LHS
            Sₖ = FractalIntegrals.singlelayer_operator_helmholtz(Γ, k)
            Sₖₕ =FractalIntegrals.discretise(Sₖ, h_mesh = 0.05, h_quad = 0.01)

            # solve discrete problem
            ϕ = [Sₖₕ \ fθ[n] for n in eachindex(recipangles)]

            # initialise matrix of ffps evaluated at different incident angles
            recip_matrix = zeros(ComplexF64, N, N)

            # construct ffps and evaluate to fill matrix
            for n in eachindex(recipangles)
                ffp = FractalIntegrals.farfield_pattern(ϕ[n], k, h_quad = 0.01)
                recip_matrix[n, :] .= ffp.(recipangles)
            end
            @testset "$Γname, wavenumber $k" begin
                @test recip_matrix ≈ transpose(recip_matrix) rtol = 1e-2
            end
        end
    end
end

@testset "Near- and far-field Comparison" begin
    θ = 0:0.01:2π
    d = [1, 1]/sqrt(2)
    r = 10 # some big number, to approximate the limit r → ∞
    x = [r*[cos(θ_), sin(θ_)] for θ_ in θ]
    k = 5
    f₁(x) = exp(im*k*d[1]*x)
    f₂(x) = exp(im*k*(d[1]*x[1]+d[2]*x[2]))

    for Γname in ["cantor set", "cantor dust"]
        Γ = getfractal(Γname)
        # data depends on ambient dimension
        Γ.n == 1 ? f = f₁ : f = f₂
        # operators, and the discretisation
        Sₖ = FractalIntegrals.singlelayer_operator_helmholtz(Γ, k)
        Sₖₕ = FractalIntegrals.discretise(Sₖ, h_mesh = 0.05, h_quad = 0.01)
        ϕ = Sₖₕ \ f 
        # get near- and far-field operators
        ffp = FractalIntegrals.farfield_pattern(ϕ, k)
        slp = FractalIntegrals.singlelayer_potential_helmholtz(ϕ, k)
        @testset "$Γname" begin
            @test (exp(im*k*r)/sqrt(r)) * ffp.(θ) ≈ slp.(x)
        end
    end
end

# other tests to implement

# quadrature
    # check gauss and barycentre converge to the same thing, for invaiant measures 