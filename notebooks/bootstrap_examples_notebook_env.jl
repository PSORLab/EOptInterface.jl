import Pkg

"""
    find_eoi_repo_root(start_dir=pwd())

Walk upward from `start_dir` until the `EOptInterface.jl` repo root is found.
The root is identified by the presence of:

- `Project.toml`
- `src/EOptInterface.jl`
- `examples/Project.toml`
- `notebooks/eoi_publication_plots.jl`
"""
function find_eoi_repo_root(start_dir::AbstractString=pwd())
    dir = abspath(start_dir)
    while true
        if isfile(joinpath(dir, "Project.toml")) &&
           isfile(joinpath(dir, "src", "EOptInterface.jl")) &&
           isfile(joinpath(dir, "examples", "Project.toml")) &&
           isfile(joinpath(dir, "notebooks", "eoi_publication_plots.jl"))
            return dir
        end
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    error(
        "Could not locate the EOptInterface.jl repo root from:\n  $(abspath(start_dir))\n" *
        "Expected Project.toml, src/EOptInterface.jl, examples/Project.toml, and notebooks/eoi_publication_plots.jl."
    )
end

"""
    prepare_examples_notebook_env(start_dir=pwd(); instantiate=true, io=devnull)

Prepare the shared `examples/` environment used by the NDMC notebook:

- locate the repo root,
- activate `examples/Project.toml`,
- instantiate missing dependencies,
- return commonly used repo paths for later notebook cells.
"""
function prepare_examples_notebook_env(
    start_dir::AbstractString=pwd();
    instantiate::Bool=true,
    io=devnull,
)
    repo_root = find_eoi_repo_root(start_dir)
    examples_env = joinpath(repo_root, "examples")

    if Base.active_project() != joinpath(examples_env, "Project.toml")
        Pkg.activate(examples_env; io=io)
    end
    "@" in LOAD_PATH || pushfirst!(LOAD_PATH, "@")
    instantiate && Pkg.instantiate(; io=io)

    return (
        repo_root = repo_root,
        examples_env = examples_env,
        generated_dir = joinpath(repo_root, "examples", "generated"),
        ndmc_module_path = joinpath(repo_root, "examples", "NDMCExample.jl"),
        style_path = joinpath(repo_root, "notebooks", "eoi_publication_plots.jl"),
    )
end
