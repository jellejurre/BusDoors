using MPI
using Gridap
using GridapGmsh
model = GmshDiscreteModel("./geometry.msh")
order = 2
u((x,y)) = (x+y)^order
f(x) = -Δ(u,x)
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags="Hand")
U = TrialFESpace(u,V)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
l(v) = ∫( v*f )dΩ
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
writevtk(Ω,"results",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
