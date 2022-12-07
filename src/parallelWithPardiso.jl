using Gridap
using Pardiso
using SparseMatricesCSR

model = GmshDiscreteModel("./geometry.msh")
order=2
u((x,y)) = (x+y)^order
f(x) = -Δ(u,x)
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,dirichlet_tags="DirichletEdges")
U = TrialFESpace(u,V)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
l(v) = ∫( v*f )dΩ
assem = SparseMatrixAssembler(SparseMatrixCSR{1,Float64,Int},Vector{Float64},U,V)
op = AffineFEOperator(a,l,U,V,assem)
ls = PardisoSolver()
solver = LinearFESolver(ls)
uh = solve(solver,op)
writevtk(Ω,"results",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])