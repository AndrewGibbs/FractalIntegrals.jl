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
							farfield_pattern,
							singlelayer_potential_helmholtz

# ╔═╡ 63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
Γ = getfractal("cantorset")

# ╔═╡ 85d36bae-021e-4a52-b767-d65e1f3fc981
S = singlelayer_operator_helmholtz(Γ, 2, ambient_dimension = 2)

# ╔═╡ e5c313e7-ae97-40de-b666-e5e854e1eea0
k = 10

# ╔═╡ 92fb4162-715a-4a4c-b0b0-9b0ee0a8146d
d = [1, 1]/sqrt(2)

# ╔═╡ 1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
f(x) = exp(im*k*d[1])

# ╔═╡ 76924a03-6acb-4e39-a826-078e684cd9da
ϕ = S \ f

# ╔═╡ bd4de407-a217-4d6a-99e0-b98f35efca37
θ = 0:0.01:2π

# ╔═╡ c6ffa209-db0f-4ff5-84e9-84ad165fc3db
ffp = farfield_pattern(ϕ, k)

# ╔═╡ 8f26187a-8dbe-4e3d-8b9b-322d2296a8da
plot(θ, [real.(ffp.(θ)), imag.(ffp.(θ))], labels = ["realffp" "imagffp"])

# ╔═╡ Cell order:
# ╠═38d67822-1d98-11ef-1aea-915089b3d1a9
# ╠═1115792d-666a-4ca6-9a0a-1a9c6df6b521
# ╠═66ab830e-b2f8-479a-b448-691d549d6f13
# ╠═63fb2b9b-1fbb-46d1-bf6f-4dafbd263871
# ╠═85d36bae-021e-4a52-b767-d65e1f3fc981
# ╠═e5c313e7-ae97-40de-b666-e5e854e1eea0
# ╟─92fb4162-715a-4a4c-b0b0-9b0ee0a8146d
# ╠═1640b341-3d3f-4a5d-b2c0-d62114ffbf4b
# ╠═76924a03-6acb-4e39-a826-078e684cd9da
# ╠═bd4de407-a217-4d6a-99e0-b98f35efca37
# ╠═c6ffa209-db0f-4ff5-84e9-84ad165fc3db
# ╠═8f26187a-8dbe-4e3d-8b9b-322d2296a8da
