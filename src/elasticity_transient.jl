# import libraries

using Gridap
using GridapGmsh
using LinearAlgebra
using Statistics: mean
using Gridap.Geometry
using Base

# Setup model
model = GmshDiscreteModel("geometry.msh")
labels = get_face_labeling(model)
dimension = 3
order = 1

# Setup fields
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)

V0 = TestFESpace(model,reffe;conformity=:H1,
    dirichlet_tags=["Hinge ceiling top", "Hinge ceiling bottom", "Area top"],
    dirichlet_masks=[(true, true, true), (true, true, true), (true,true,true)])

g1(x, t::Real) = VectorValue(0,0,0)
g1(t::Real) = x -> g1(x,t)

U = TransientTrialFESpace(V0, [g1,g1,g1])

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

function normalzed_bell_curve(mean, deviation, x)
    return  ℯ^(-norm(x - mean)/(2 * deviation^2))
end

f(t, x) = normalzed_bell_curve(5, 1, t) * amplitude*exp(-deviation*norm(x-x0))*direction
f(t) = x -> f(t, x) 

# Set up equations for Transient case

degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

m(utt,v) = ∫(v⊙utt)dΩ

c(ut,v) = ∫(v⊙ut)dΩ

a(u,v) = ∫(ε(v) ⊙ (σ_bimat∘(ε(u),tags)))dΩ 

b(t,v) = ∫(v⋅f(t))dΩ
    
m(t,utt,v) = m(utt,v)
c(t,ut,v) = c(ut,v)
a(t,u,v) = a(u,v)

res(t,u,v) = m(∂tt(u),v) + c(∂t(u),v) + a(u,v) - b(t,v)
jac(t,u,du,v) = a(du,v)
jac_t(t,u,dut,v) = c(dut,v)
jac_tt(t,u,dutt,v) = m(dutt,v)


op = TransientFEOperator(res, jac, jac_t, jac_tt, U, V0)
# op = TransientFEOperator(res, U, V0; order=2) # for auto differentiation, without the need to provide jacobians

linear_solver = LUSolver()

dt = 0.5

# with γ = 0.5 and β = 0.25 the newmark method is unconditionally stable with any dt

γ = 0.5
β = 0.25

u₀ = interpolate_everywhere(VectorValue(0,0,0),U(0))
v₀ = interpolate_everywhere(VectorValue(0,0,0), U(0))
a₀ = interpolate_everywhere(VectorValue(0,0,0), U(0))

ode_solver = Newmark(linear_solver, dt, γ, β)

t₀ = 0.0
T = 10.0

uₕₜ = solve(ode_solver,op,(u₀,v₀, a₀),t₀,T)

createpvd("transient_elasticity_results") do pvd
    pvd[0] = createvtk(Ω,"result_transient0.0.vtu",cellfields=["u"=>u₀, "f"=>f(0), 
    "vonmises"=>σ_vm∘(σ_bimat∘(ε(u₀),tags))])
    for (uₕ,t) in uₕₜ
        pvd[t] = createvtk(Ω,"result_transient$t"*".vtu",cellfields=["u"=>uₕ, "f"=>f(t), 
                "vonmises"=>σ_vm∘(σ_bimat∘(ε(uₕ),tags))])
    end
end