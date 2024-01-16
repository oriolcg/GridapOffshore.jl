module Timoshenko
using Gridap

# Discrete model
L = 3.0
domain = (0,L)
partition = (10,)
model = CartesianDiscreteModel(domain,partition)

# Boundary labels
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"left",[1])
add_tag_from_tags!(labels,"right",[2])

# Triangulations
Ω = Triangulation(model)
Γₗ = Boundary(Ω,tags="left")
Γᵣ = Boundary(Ω,tags="right")

# FE spaces
order = 1
reffeᵥ = ReferenceFE(lagrangian,Float64,order)
reffeθ = ReferenceFE(lagrangian,VectorValue{1,Float64},order)
Vₕ = TestFESpace(Ω,reffeᵥ,dirichlet_tags="left")
Ψₕ = TestFESpace(Ω,reffeθ,dirichlet_tags="left")
Uₕ = TrialFESpace(Vₕ,0.0)
Θₕ = TrialFESpace(Ψₕ,VectorValue(0.0))
Xₕ = MultiFieldFESpace([Uₕ,Θₕ])
Yₕ = MultiFieldFESpace([Vₕ,Ψₕ])

# Integration measures and normals
dΩ = Measure(Ω,2*order)
dΓᵣ = Measure(Γᵣ,2*order)

# Problem parameters
b =1.0
E = 1.0e6
ν = 0.2
G = E/(2*(1+ν))
EI(t) = b/12*t^3*E
GA(t) = G*t*b
F = -1.0

# Weak form
a(t,(v,θ),(w,ψ)) = ∫( EI(t)*(∇(θ)⊙∇(ψ)) + GA(t)*((∇(v)-θ)⋅(∇(w)-ψ)) )dΩ
l((w,ψ)) = ∫( F*w )dΓᵣ
op(t) = AffineFEOperator((x,y)->a(t,x,y),l,Xₕ,Yₕ)

# Solution
t₀ = 0.1
vₕ,θₕ = solve(op(t₀))

# Postprocessing
writevtk(Ω,"Timoshenko_solution",cellfields=["v"=>vₕ,"theta"=>θₕ])
xₚ = Point((L,))
vₑ(t) = F*L^3/(3*EI(t))
println("Computed solution: $(vₕ(xₚ)), Exact solution: $(vₑ(t₀)), Error: $(abs(vₕ(xₚ)-vₑ(t₀)))")

# Shear Locking
ts = Float64[]
v_rel = Float64[]
for i in 0:10
  tᵢ = 2.0^(-float(i))
  vₕ,θₕ = solve(op(tᵢ))
  push!(ts,tᵢ)
  push!(v_rel,abs(vₕ(xₚ)/vₑ(tᵢ)))
end

# Plot results
using Plots
plt = plot(xlabel="t",ylabel="vₕ/vₑ",title="Shear Locking",legend=:bottomright,xscale=:log10,yscale=:log10)
plot!(plt,ts,v_rel,marker=:circle,label="Galerkin")
plot!(plt,ts,ones(length(ts)),ls=:dash,color=:black,label=false)
display(plt)

# Solution 1: subintegration
dΩ₁ = Measure(Ω,1)
a₁(t,(v,θ),(w,ψ)) = ∫( GA(t)*(∇(v)⋅∇(w)) )dΩ + ∫( EI(t)*(∇(θ)⊙∇(ψ)) + GA(t)*((-∇(v)⋅ψ-θ⋅∇(w)+θ⋅ψ)) )dΩ₁
op₁(t) = AffineFEOperator((x,y)->a₁(t,x,y),l,Xₕ,Yₕ)
v_rel₁ = Float64[]
for tᵢ in ts
  vₕ,θₕ = solve(op₁(tᵢ))
  push!(v_rel₁,abs(vₕ(xₚ)/vₑ(tᵢ)))
end

# Plot results with subintegration
plot!(plt,ts,v_rel₁,marker=:square,label="Sub-integration")
display(plt)

# Solution 2: compatible spaces
order = 2
reffeᵥ₂ = ReferenceFE(lagrangian,Float64,order)
reffeθ₂ = ReferenceFE(lagrangian,VectorValue{1,Float64},order-1)
Vₕ₂ = TestFESpace(Ω,reffeᵥ₂,dirichlet_tags="left")
Ψₕ₂ = TestFESpace(Ω,reffeθ₂,dirichlet_tags="left")
Uₕ₂ = TrialFESpace(Vₕ₂,0.0)
Θₕ₂ = TrialFESpace(Ψₕ₂,VectorValue(0.0))
Xₕ₂ = MultiFieldFESpace([Uₕ₂,Θₕ₂])
Yₕ₂ = MultiFieldFESpace([Vₕ₂,Ψₕ₂])
op₂(t) = AffineFEOperator((x,y)->a(t,x,y),l,Xₕ₂,Yₕ₂)
v_rel₂ = Float64[]
for tᵢ in ts
  vₕ,θₕ = solve(op₂(tᵢ))
  push!(v_rel₂,abs(vₕ(xₚ)/vₑ(tᵢ)))
end

# Plot results with compatible spaces
plot!(plt,ts,v_rel₂,marker=:diamond,label="Q2/Q1")
display(plt)

# Solution 3: Stabilized formulation
h = L/partition[1]
c₁ = 12.0; c₂ = 1.0; c₃ = 12.0; c₄ = 1.0
τᵥ(t) = 1.0/(c₃*GA(t)/(h^2)+c₄*(GA(t)^2/EI(t)))
τθ(t) = 1.0/(c₁*EI(t)/(h^2)+c₂*GA(t))
a₃(t,(v,θ),(w,ψ)) = ∫( (EI(t)-τᵥ(t)*GA(t)^2)*(∇(θ)⊙∇(ψ)) + GA(t)*(1-τθ(t)*GA(t))*((∇(v)-θ)⋅(∇(w)-ψ)) )dΩ
op₃(t) = AffineFEOperator((x,y)->a₃(t,x,y),l,Xₕ,Yₕ)
v_rel₃ = Float64[]
for tᵢ in ts
  vₕ,θₕ = solve(op₃(tᵢ))
  push!(v_rel₃,abs(vₕ(xₚ)/vₑ(tᵢ)))
end

# Plot results with stabilized formulation
plot!(plt,ts,v_rel₃,marker=:utriangle,label="Stabilized")
display(plt)
end
