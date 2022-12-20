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

function main(parts)
  options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  GridapPETSc.with(args=split(options)) do
    model = GmshDiscreteModel("./geometry.msh")
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
    # LU Solver
    time = @elapsed begin 
      uh = solve(op)
    end
    @printf("LU Solver: %.20f", time)
    # Backslash Solver
    op = AffineFEOperator(a,l,U,V)
    solver = LinearFESolver(BackslashSolver())
    time = @elapsed begin 
      uh = solve(solver, op)
    end
    @printf("Backslash Solver: %.20f", time)
    # Newton Raphson Solver
    op = AffineFEOperator(a,l,U,V)
    solver = FESolver()
    time = @elapsed begin 
      uh = solve(solver, op)
    end
    @printf("Newton Raphson Solver: %.20f", time)
    # PETScLinearSolver
    solver = PETScLinearSolver()
    time = @elapsed begin 
      uh = solve(solver, op)
    end
    @printf("PETScLinearSolver: %.20f", time)
    # PETScNonlinearSolver
    # solver = PETScNonlinearSolver()
    # time = @elapsed begin 
    #   uh = solve(solver, op)
    # end
    @printf("PETScNonlinearSolver: %.20f", time)
    writevtk(Ω,"results",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
  end
end

# Run the main funcion in a parallel fashion
prun(main, mpi, n)
