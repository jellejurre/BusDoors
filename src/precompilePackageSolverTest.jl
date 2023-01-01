using PackageCompiler
PackageCompiler.create_sysimage(
    ["MPI", "Gridap", "GridapGmsh", "GridapPETSc", "GridapDistributed", "PartitionedArrays", "Printf"]; 
    sysimage_path="PreCompileSolverTest.so", precompile_execution_file="./src/poisson.jl")