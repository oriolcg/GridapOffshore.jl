module Euler_Bernoulli
using Gridap

# Discrete model
L = 3.0
domain = (0,L)
partition = (100,)
model = CartesianDiscreteModel(domain,partition)

# Boundary labels
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"left",[1])
add_tag_from_tags!(labels,"right",[2])

# Triangulations
Ω = Triangulation(model)
Γₗ = Boundary(Ω,tags="left")
Γᵣ = Boundary(Ω,tags="right")
Λ = Skeleton(Ω)

# Output triangulations to VTK files
writevtk(Ω,"Omega")
writevtk(Γₗ,"Gamma_left")
writevtk(Γᵣ,"Gamma_right")
writevtk(Λ,"Lambda")

# FE spaces
order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
Vₕ = TestFESpace(Ω,reffe,dirichlet_tags="left")
Uₕ = TrialFESpace(Vₕ,0.0)

# Integration measures and normals
dΩ = Measure(Ω,2*order)
dΓₗ = Measure(Γₗ,1)
dΛ = Measure(Λ,1)
nΓₗ = get_normal_vector(Γₗ)
nΛ = get_normal_vector(Λ)

# Problem parameters
EI = 1.5e7
q(x) = -1.0e3*(L-x[1]) # triangular dustributed load

# Weak form
h = L/partition[1]
a(v,w) = ∫( EI*(Δ(v)*Δ(w)) )dΩ +
  ∫( 2*EI/h*((∇(v)⋅nΓₗ)*(∇(w)⋅nΓₗ)) -
    (EI*Δ(v))*(∇(w)⋅nΓₗ) -
    (EI*Δ(w))*(∇(v)⋅nΓₗ) )dΓₗ +
  ∫( EI/h*(jump(∇(v)⋅nΛ)*jump(∇(w)⋅nΛ)) -
    mean(EI*Δ(v))*jump(∇(w)⋅nΛ) -
    mean(EI*Δ(w))*jump(∇(v)⋅nΛ) )dΛ
l(w) = ∫( q*w )dΩ
op = AffineFEOperator(a,l,Uₕ,Vₕ)

# Solution
vₕ = solve(op)

#Postprocessing
writevtk(Ω,"EB_solution",cellfields=["v"=>vₕ,"theta"=>∇(vₕ)])
xₚ = Point((L,))
println(vₕ(xₚ))

end
