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
