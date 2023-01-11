using PackageCompiler
PackageCompiler.create_sysimage(
    ["Gridap", "GridapDistributed", "GridapPETSc", "GridapGmsh", "MPI", "PartitionedArrays"]; 
    sysimage_path="PreCompileGenerate.so", precompile_execution_file="./src/generateMesh.jl")