using Documenter
using DocumenterTools
using EOptInterface

makedocs(;
    modules = [EOptInterface],
    authors = "Joseph Choi, Dimitri Alston, Pengfei Xu, and Matthew Stuber",
    sitename = "EOptInterface.jl",
    format = Documenter.HTML(;
             canonical = "https://PSORLab.github.io/EOptInterface.jl/stable",
             collapselevel = 1,
             assets = ["assets/favicon.ico"],
    ),
    pages = Any[
        "API Reference" => "index.md",
        "News" => "news.md",
    ],
)

deploydocs(;
    repo = "github.com/PSORLab/EOptInterface.jl",
)
