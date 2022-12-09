using Gridap
using GridapGmsh
import LinearAlgebra: norm
model = GmshDiscreteModel("./geometry.msh")

# order = 1

# reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
# V0 = TestFESpace(model,reffe;
#   conformity=:H1,
#   dirichlet_tags="DirichletEdges",
#   dirichlet_masks=[(true,false,false), (true,true,true)])
# g1(x) = VectorValue(0.005,0.0,0.0)
# g2(x) = VectorValue(0.0,0.0,0.0)
# U = TrialFESpace(V0,[g1,g2])
# const E = 70.0e9
# const ν = 0.33
# const λ = (E*ν)/((1+ν)*(1-2*ν))
# const μ = E/(2*(1+ν))
# σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
# degree = 2*order
# Ω = Triangulation(model)
# dΩ = Measure(Ω,degree)
# a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
# l(v) = 0
# op = AffineFEOperator(a,l,U,V0)
# uh = solve(op)
# writevtk(Ω,"results_elasticity",cellfields=["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ∘ε(uh)])

using Gridap.Geometry
labels = get_face_labeling(model)
dimension = 3
tags = get_face_tag(labels,dimension)
const glass_tag = get_tag_from_name(labels,"Glass")
const aluminium_tag = get_tag_from_name(labels,"Aluminium")

order = 1

reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
# V0 = TestFESpace(model,reffe;
#   conformity=:H1,
#   dirichlet_tags="DirichletEdges",
#   dirichlet_masks=[(true,false,false), (true,true,true)]
#   )
V0 = TestFESpace(model,reffe;
  conformity=:H1
  )
# g1(x) = VectorValue(0.005,0.0,0.0)
# g2(x) = VectorValue(0.0,0.0,0.0)
# U = TrialFESpace(V0,[g1,g2])
U = TrialFESpace(V0)

function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

const E_alu = 70.0e9  # Young's modulus
const ν_alu = 0.33  # poisson ratio
const (λ_alu,μ_alu) = lame_parameters(E_alu,ν_alu)

const E_glass = 2855
const ν_glass = 0.35
const (λ_glass,μ_glass) = lame_parameters(E_glass,ν_glass)

function σ_bimat(ε,tag)
  if tag == glass_tag
    return λ_glass*tr(ε)*one(ε) + 2*μ_glass*ε
  else
    return λ_alu*tr(ε)*one(ε) + 2*μ_alu*ε
  end
end

sources = [(10, 1, VectorValue(0.6,1,0))]

function get_force(x, tag)
  # if (tag != aluminium_tag)
  #   return 0
  # end
  sum_forces = 0
  for source in sources
    val = source[1] * ℯ^(-(0.1 + source[2] * norm(source[3] - x)))
    sum_forces += val
  end
  return sum_forces
end

σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)

a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags)) )*dΩ
l(v) = ∫( get_force∘(v, tags) )*dΓ

op = AffineFEOperator(a,l,U,V0)
uh = solve(op)
writevtk(Ω,"results_bimat",cellfields=
  ["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ_bimat∘(ε(uh),tags)])
