using PackageCompiler
PackageCompiler.create_sysimage(
    ["Gridap", "GridapDistributed", "GridapPETSc", "GridapGmsh", "MPI", "PartitionedArrays"]; 
    sysimage_path="PreCompileRun.so", precompile_execution_file="./src/poisson.jl")