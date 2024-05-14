using FractalIntegrals
using Test

@testset "Preset fractals" begin
    for Γ in keys(FractalIntegrals.fractaldict)
        strΓ = String(Γ)
        @testset "$strΓ" begin
            @test_nowarn FractalIntegrals.getfractal(Γ)
            @test_nowarn FractalIntegrals.getfractal(BigFloat, Γ)
        end
    end
end

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

# other tests to implement

# wave-based
    # reciprocity
    # check that scattered field converges to far-field
    # test against Andrea's old code

# quadrature
    # check s-energy line integrals for singularities which can be calculated by hand
    # check gauss and barycentre converge to the same thing, for invaiant measures 