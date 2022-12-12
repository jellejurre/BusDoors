using Gridap
using GridapGmsh
import LinearAlgebra: norm
import Statistics: mean
using Gridap.Geometry

# Setup model
model = GmshDiscreteModel("./geometry.msh")
labels = get_face_labeling(model)
dimension = 3
order = 1

# Setup fields
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V0 = TestFESpace(model,reffe;conformity=:H1)
U = TrialFESpace(V0)

# Setup tags
tags = get_face_tag(labels,dimension)
const glass_tag = get_tag_from_name(labels,"Glass")
const aluminium_tag = get_tag_from_name(labels,"Aluminium")

# Setup lamé parameters
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

# Setup forces
sources = [(10, 1, VectorValue(0,0,0), VectorValue(0, 0, 1))]

function bell_curve(mean, deviation, x)
  return 1/(deviation * sqrt(2 * π)) * ℯ^(-norm(x - mean)/(2 * deviation^2))
end

function get_force(x)
  nodeAverage = mean(x) #Average required because we have the coordinates of the four corners of the cell
  sum_forces = VectorValue(0, 0, 0)
  for source in sources
    amplitude = source[1] * bell_curve(source[3], source[2], nodeAverage)
    force_val = amplitude * source[4]
    sum_forces += force_val
  end
  return sum_forces
end

forces = map((cell) -> get_force(cell), get_cell_coordinates(model))

# Set up equation
σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags)) )*dΩ
l(v) = ∫( v⋅forces )*dΩ

# Perform equation
op = AffineFEOperator(a,l,U,V0)
uh = solve(op)

# Write output
writevtk(Ω,"results",cellfields=
  ["force"=>forces, "uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ_bimat∘(ε(uh),tags)])