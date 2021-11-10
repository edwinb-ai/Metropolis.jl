using Metropolis
using Documenter

makedocs(;
    modules=[Metropolis],
    authors="Edwin Bedolla",
    repo="https://github.com/edwinb-ai/Metropolis.jl/blob/{commit}{path}#L{line}",
    sitename="Metropolis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://edwinb-ai.github.io/Metropolis.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/edwinb-ai/Metropolis.jl")
