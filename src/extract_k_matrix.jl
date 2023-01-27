using Gridap
using GridapGmsh
using LinearAlgebra
using Statistics: mean
using Gridap.Geometry
using Base
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Fields
using Gridap.CellData
using FillArrays
using Test
using InteractiveUtils

model = GmshDiscreteModel("geometry.msh")
Th = Triangulation(model)
order = 1
dimension = 3
Qh = CellQuadrature(Th,4*order)
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

V0 = TestFESpace(model,reffe;conformity=:H1,
    dirichlet_tags=["Hinge ceiling top", "Hinge ceiling bottom", "Area top"],
    dirichlet_masks=[(true, true, true), (true, true, true), (true,true,true)])
g1(x) = VectorValue(0,0,0)
U0 = TrialFESpace(V0, [g1,g1,g1])

labels = get_face_labeling(model)
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

du = get_trial_fe_basis(U0)
dv = get_fe_basis(V0)
jac = integrate(ε(dv) ⊙ (σ_bimat∘(ε(du),tags)),Qh)
sigmak = get_cell_dof_ids(U0)
assem = SparseMatrixAssembler(U0,V0);
jac = collect(jac)
rs = ([jac],[sigmak],[sigmak])
A = allocate_matrix(assem,rs)
A = assemble_matrix!(A,assem,rs)

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