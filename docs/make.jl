using SymPhase
using Documenter

DocMeta.setdocmeta!(SymPhase, :DocTestSetup, :(using SymPhase); recursive=true)

makedocs(;
    modules=[SymPhase],
    authors="njuwfang <njuwfang@gmail.com> and contributors",
    repo="https://github.com/njuwfang/SymPhase.jl/blob/{commit}{path}#{line}",
    sitename="SymPhase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://njuwfang.github.io/SymPhase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/njuwfang/SymPhase.jl",
    devbranch="main",
)
