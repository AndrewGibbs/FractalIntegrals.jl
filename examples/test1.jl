### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 38d67822-1d98-11ef-1aea-915089b3d1a9
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
end

# ╔═╡ 1115792d-666a-4ca6-9a0a-1a9c6df6b521
using FractalIntegrals, Plots

# ╔═╡ 66ab830e-b2f8-479a-b448-691d549d6f13
import FractalIntegrals: 	singlelayer_operator_laplace,
							HausdorffMeasure,
							singlelayer_operator_helmholtz,
							farfield_pattern

# ╔═╡ 63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
Γ = getfractal("cantorset")

# ╔═╡ d5dc2dd8-565d-45cd-9467-b1724761a1a1
plot(Γ)

# ╔═╡ 23858b6e-db90-4249-a71c-0fe28dc1bdc4
ℋ = HausdorffMeasure(Γ)

# ╔═╡ 89d50860-c0d9-4bdb-b014-b8f7a9b6ed6d
Slap = singlelayer_operator_helmholtz(Γ, ambient_dimension = 2)

# ╔═╡ 85d36bae-021e-4a52-b767-d65e1f3fc981
S = singlelayer_operator_helmholtz(ℋ, 2, ambient_dimension = 2)

# ╔═╡ 1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
f(x) = exp(im*5*x)

# ╔═╡ 76924a03-6acb-4e39-a826-078e684cd9da
ϕ = S \ f

# ╔═╡ 1c5da493-2e90-4731-91ef-e6fe4aae1031
ffp = farfield_pattern(ϕ, 5)

# ╔═╡ Cell order:
# ╠═38d67822-1d98-11ef-1aea-915089b3d1a9
# ╠═1115792d-666a-4ca6-9a0a-1a9c6df6b521
# ╠═66ab830e-b2f8-479a-b448-691d549d6f13
# ╠═63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
# ╠═d5dc2dd8-565d-45cd-9467-b1724761a1a1
# ╠═23858b6e-db90-4249-a71c-0fe28dc1bdc4
# ╠═89d50860-c0d9-4bdb-b014-b8f7a9b6ed6d
# ╠═85d36bae-021e-4a52-b767-d65e1f3fc981
# ╠═1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
# ╠═76924a03-6acb-4e39-a826-078e684cd9da
# ╠═1c5da493-2e90-4731-91ef-e6fe4aae1031
