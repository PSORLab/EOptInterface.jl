using EOptInterface
using Documenter

DocMeta.setdocmeta!(EOptInterface, :DocTestSetup, :(using EOptInterface); recursive=true)

makedocs(;
    modules=[EOptInterface],
    authors="Joseph Choi, Dimitri Alston, Pengfei Xu, and Matthew Stuber",
    sitename="EOptInterface.jl",
    format=Documenter.HTML(;
        canonical="https://PSORLab.github.io/EOptInterface.jl/mpc/",
        edit_link="MPC",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "MPC Module" => "mpc_module.md",
    ],
)

deploydocs(;
    repo="github.com/PSORLab/EOptInterface.jl",
    devbranch="MPC",
    devurl="mpc",
)
