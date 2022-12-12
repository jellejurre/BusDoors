using MPI
using Gridap
using GridapGmsh
using GridapPETSc
using GridapDistributed
using PartitionedArrays

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
    uh = solve(op)
    writevtk(Ω,"results",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
  end
end

# Run the main funcion in a parallel fashion
prun(main, mpi, n)
