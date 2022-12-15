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
sources = [(10, 5, VectorValue(0.3,1,0), VectorValue(0, 0, 1))]
force_tag = get_tag_from_name(labels,"Glass")

function bell_curve(mean, deviation, x)
  return 1/(deviation * sqrt(2 * π)) * ℯ^(-norm(x - mean)/(2 * deviation^2))
end


function get_volume(points)
  point1 = points[1]
  point2 = points[2]
  point3 = points[3]
  point4 = points[4]
  vec1to2 = point2-point1
  vec1to3 = point3-point1
  vec1to4 = point4-point1
  return (1/6) * abs(dot(cross(vec1to2, vec1to3), vec1to4)) 
end

force_tag_cells = Triangulation(model, tags=[force_tag])
total_force_volume = sum(map((cell) -> get_volume(cell), get_cell_coordinates(force_tag_cells)))

function get_force(x)
  volume_fraction = get_volume(x)/total_force_volume
  nodeAverage = mean(x) #Average required because we have the coordinates of the four corners of the cell
  sum_forces = VectorValue(0, 0, 0)
  for source in sources
    amplitude = source[1] * bell_curve(source[3], source[2], nodeAverage)
    force_val = amplitude * source[4]
    normalized_force_val = volume_fraction * force_val
    sum_forces += normalized_force_val
  end
  return sum_forces
end

forces = map((cell) -> get_force(cell), get_cell_coordinates(model))

function filter_tags(v, force, tag)
  if tag == force_tag
    return force ⋅ v
  else 
    return 0.0
  end
end

# Set up equation
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags)) )*dΩ
l(v) = ∫( filter_tags∘(v,forces, tags) )*dΩ

# Perform equation
op = AffineFEOperator(a,l,U,V0)
uh = solve(op)

# Write output
writevtk(Ω,"results",cellfields=
  ["force"=>forces, "uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ_bimat∘(ε(uh),tags)])