# BusDoors

## How to run
You can run the julia files with the `julia --project=. [file to run]` command. This will use the Project.toml and Manifest.toml to do dependency resolution. If you want to add new dependencies to the project, you can do so by typing `activate .` in the Julia REPL package manager, and any new packages you add will be added to the .toml files.

If you don't have the packages installed yet, you can install the required packages by typing the `activate .` and `instantiate` commands in the Julia REPL package manager.