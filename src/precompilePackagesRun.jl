using PackageCompiler
PackageCompiler.create_sysimage(
    ["Gridap", "GridapGmsh"]; 
    sysimage_path="PreCompileRun.so", precompile_execution_file="./src/elasticity_transient.jl")