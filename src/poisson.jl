using MPI
using Gridap
using GridapGmsh
model = GmshDiscreteModel("./geometry_beam.msh")
order = 2
u((x,y)) = (x+y)^order
f(x) = -Δ(u,x)
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags="Dirichlet1")
U = TrialFESpace(u,V)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
l(v) = ∫( v*f )dΩ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
writevtk(Ω,"results_beam",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
