using MPI
using Gridap
using GridapGmsh
using GridapPETSc
using GridapDistributed
using PartitionedArrays
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
n = MPI.Comm_size(comm)
rank= MPI.Comm_rank(comm)

function main(parts)
  options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  GridapPETSc.with(args=split(options)) do
    read_model_time = @elapsed begin
      model = GmshDiscreteModel(parts, "./geometry.msh")
    end

    precompute_time = @elapsed begin
      order = 2
      u((x,y)) = (x+y)^order
      f(x) = -Δ(u,x)
      reffe = ReferenceFE(lagrangian,Float64,order)
      V = TestFESpace(model,reffe,dirichlet_tags="DirichletEdges")
      U = TrialFESpace(u,V)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,2*order)
      a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
      l(v) = ∫( v*f )dΩ
      op = AffineFEOperator(a,l,U,V)
      backslash_solver = BackslashSolver()
      PETSC_solver = PETScLinearSolver()
    end

    LU_compute_time = @elapsed begin
      uh = @time solve(op)
    end

    Backslash_compute_time = @elapsed begin
      uh = @time solve(backslash_solver, op)
    end

    PETSC_compute_time = @elapsed begin
      uh = @time solve(PETSC_solver, op)
    end

    if (rank == 0)
      @printf "Time for reading model: %.5f\n" read_model_time
      @printf "Time for precomputation: %.5f\n" precompute_time
      @printf "Time for LU Solver: %.5f\n" LU_compute_time
      @printf "Time for Backslash Solver: %.5f\n" Backslash_compute_time
      @printf "Time for PETSC Solver: %.5f\n" PETSC_compute_time
    end
    writevtk(Ω,"results",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
  end
end

# Run the main funcion in a parallel fashion
prun(main, mpi, n)
