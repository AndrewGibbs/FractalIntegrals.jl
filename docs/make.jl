push!(LOAD_PATH,"../src/")
using FractalIntegrals
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
         sitename = "FractalIntegrals.jl",
         modules  = [FractalIntegrals],
         checkdocs=:none,
         pages=[
                "Home" => "index.md",
                "Fractals" => "makeIFS.md",
                "Integrals" => "integrals.md",
                "References" => "refs.md"
               ],
        format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
        plugins=[bib]
        )

deploydocs(;
    repo="github.com/AndrewGibbs/FractalIntegrals.jl",
)