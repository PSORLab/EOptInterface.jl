using EOptInterface
using Documenter

DocMeta.setdocmeta!(EOptInterface, :DocTestSetup, :(using EOptInterface); recursive=true)

makedocs(;
    modules=[EOptInterface],
    authors="Joseph Choi <joseph03choi@gmail.com>",
    sitename="EOptInterface.jl",
    format=Documenter.HTML(;
        canonical="https://joseph03choi.github.io/EOptInterface.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/joseph03choi/EOptInterface.jl",
    devbranch="master",
)
