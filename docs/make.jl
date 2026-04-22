using Documenter
using MarineShelf
using OceanInterp

makedocs(
    modules  = [MarineShelf, OceanInterp],
    sitename = "MarineShelf.jl",
    authors  = "Alexander Robinson and contributors",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://fesmc.github.io/MarineShelf.jl",
        edit_link  = "main",
    ),
    pages = [
        "Home"         => "index.md",
        "User Guide"   => "guide.md",
        "API Reference" => "api.md",
        "OceanInterp"  => "oceaninterp.md",
    ],
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

deploydocs(
    repo       = "github.com/fesmc/MarineShelf.jl",
    devbranch  = "main",
    push_preview = true,
)
