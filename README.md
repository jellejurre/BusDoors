# BusDoors

## How to run
You can run the julia files with the `julia --project=. [file to run]` command. This will use the Project.toml and Manifest.toml to do dependency resolution. If you want to add new dependencies to the project, you can do so by typing `activate .` in the Julia REPL package manager, and any new packages you add will be added to the .toml files.

If you don't have the packages installed yet, you can install the required packages by typing the `activate .` and `instantiate` commands in the Julia REPL package manager.

## Run in parallel
To run the FEM process in parallel, MPI need to be installed. `Microsoft MPI` for Windows, `MPICH` for other platforms.
Then use command `mpiexec -n [threadcount] julia --project=. [file to run]` to run the parallel Julia file.
`mpiexec -n [threadcount] julia --project=. [file to run]` is used for Unix-based system.

## Precompile libraries
To generate a local precompilation of the libraries for speedup, run the precompilePackageGenerate.jl and precompilePackageRun.jl files. This will create a .so file, which you can add with the `-J "[filename].so"` option. So an example final command would be `mpiexec -n [threadcount] julia -J "PreCompileRun.so" --project=. .\src\poisson.jl`