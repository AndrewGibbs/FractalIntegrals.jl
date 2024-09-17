push!(LOAD_PATH,"../src/")
using FractalIntegrals
using Documenter

makedocs(
         sitename = "FractalIntegrals.jl",
         modules  = [FractalIntegrals],
         checkdocs=:none,
         pages=[
                "Home" => "index.md",
                "Fractals" => "makeIFS.md",
                "Integrals" => "integrals.md"
               ],
               format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true"))
deploydocs(;
    repo="github.com/AndrewGibbs/FractalIntegrals.jl",
)
