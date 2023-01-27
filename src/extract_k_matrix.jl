# import libraries

using Gridap
using GridapGmsh
using LinearAlgebra
using Statistics: mean
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Algebra
using Base

# Setup model
model = GmshDiscreteModel("geometry.msh")
labels = get_face_labeling(model)
dimension = 3
order = 1

# Setup fields
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

test = TestFESpace(model,reffe;conformity=:H1,
    dirichlet_tags=["Hinge ceiling top", "Hinge ceiling bottom", "Area top"],
    dirichlet_masks=[(true, true, true), (true, true, true), (true,true,true)])

g1(x) = VectorValue(0,0,0)

trial = TrialFESpace(test, [g1,g1,g1])

# Setup tags
tags = get_face_tag(labels,dimension)
const aluminium_tag = get_tag_from_name(labels,"Aluminium")
const rubber_tag = get_tag_from_name(labels, "Rubber")
const glass_tag = get_tag_from_name(labels, "Glass")
const hinge_top_tag = get_tag_from_name(labels, "Hinge top")
const hinge_bottom_tag = get_tag_from_name(labels, "Hinge bottom")

# Setup lamé parameters, stress tensor, von mises stress and constants

function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

const E_alu = 70.0e9  # Young's modulus
const ν_alu = 0.33  # poisson ratio
const (λ_alu,μ_alu) = lame_parameters(E_alu,ν_alu)

const E_glass = 2.855e9
const ν_glass = 0.35
const (λ_glass,μ_glass) = lame_parameters(E_glass,ν_glass)

const E_rubber = 0.1e9
const ν_rubber = 0.47
const (λ_rubber, μ_rubber) = lame_parameters(E_rubber, ν_rubber)

const E_hinge = 0.01e9
const ν_hinge = 0.49
const (λ_hinge, μ_hinge) = lame_parameters(E_hinge, ν_hinge)

function σ_bimat(ε,tag)
    if tag == glass_tag
        return λ_glass*tr(ε)*one(ε) + 2*μ_glass*ε
    elseif tag == rubber_tag
        return λ_rubber*tr(ε)*one(ε) + 2*μ_rubber*ε
    elseif tag == hinge_top_tag || tag == hinge_bottom_tag 
        return λ_hinge*tr(ε)*one(ε) + 2*μ_hinge*ε
    else
        return λ_alu*tr(ε)*one(ε) + 2*μ_alu*ε
  end
end

# von misess stress
σ_vm(σ) = sqrt( 0.5 * ( (σ[1,1] - σ[2,2] )^2 + (σ[2,2]-σ[3,3])^2 + (σ[3,3]-σ[1,1])^2 ) + 
    3*(σ[1,2]^2 + σ[2,3]^2 + σ[1,3]^2) )

# setup forcing term

x0 = VectorValue(0.3, 1.3, 0.0145)
amplitude = 1000
deviation = 10 
direction = VectorValue(0,0,-1)

f(x) = amplitude*exp(-deviation*norm(x-x0))*direction

# Set up static case equations

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ε(v) ⊙ (σ_bimat∘(ε(u),tags)) )*dΩ
l(v) = ∫(f⋅v) *dΩ

# Perform equation
u = get_trial_fe_basis(trial)
v = get_fe_basis(test)
assem = SparseMatrixAssembler(trial, test)
uhd = zero(trial)
matcontribs = a(u, v)
veccontribs = l(v)
data = collect_cell_matrix_and_vector(trial,test,matcontribs,veccontribs,uhd)
A,b = assemble_matrix_and_vector(assem,data)
uh = FEFunction(trial, A\b)

for i in findall(!iszero, A)
  if abs(A[i[1],i[2]]-A[i[2],i[1]])>0.00001 
      println("Not symmetric")
  end
end

for eig in eigvals(Matrix(A))
    if eig<=0
        println("Not positive definite")
    end
end


# Write output for static case
writevtk(Ω,"static_elasticity_result",cellfields=
  ["force"=>f, "uh"=>uh,
        "epsi"=>ε(uh),
        "sigma"=>σ_bimat∘(ε(uh),tags), 
        "vonmises"=>σ_vm∘(σ_bimat∘(ε(uh),tags))])