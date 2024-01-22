module Poisson_weakBCs
using Gridap

# Define the model
L = 3.0
domain = (0.0, L, 0.0, L)
partition = (10, 10)
model = CartesianDiscreteModel(domain, partition)

# Define boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "left", [1,3,7])
add_tag_from_tags!(labels, "right", [2,4,8])
add_tag_from_tags!(labels, "bottom", [5])
add_tag_from_tags!(labels, "top", [6])

# Triangulations
Ω = Interior(model)
Γl = Boundary(model, tags="left")
Γr = Boundary(model, tags="right")
Γb = Boundary(model, tags="bottom")
Γt = Boundary(model, tags="top")

# Boundary conditions
u₀(x) = 0.0
u₁(x) = x[2]*(L-x[2])

# FE spaces
order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
Vₕ₀ = TestFESpace(model, reffe,dirichlet_tags=["left","right","bottom","top"])
Uₕᵤ = TrialFESpace(Vₕ₀,[u₁,u₁,u₀,u₀])

# Measures
degree = 2*order
dΩ = Measure(Ω, degree)

# Weak formulation strong Dirichlet
f(x) = 1.0
a(u,v) = ∫(∇(v)⋅∇(u))dΩ
b(v) = ∫(f*v)dΩ
op₀ = AffineFEOperator(a,b,Uₕᵤ,Vₕ₀)

# solve
uₕ₀ = solve(op₀)

# FE spaces with no strong Dirichlet
Vₕ = TestFESpace(model, reffe)
Uₕ = TrialFESpace(Vₕ)

# Additional measures on boundaries
dΓl = Measure(Γl, degree)
dΓr = Measure(Γr, degree)
dΓb = Measure(Γb, degree)
dΓt = Measure(Γt, degree)

# Normal vectors
nl = get_normal_vector(Γl)
nr = get_normal_vector(Γr)
nb = get_normal_vector(Γb)
nt = get_normal_vector(Γt)

# Weak formulation weak Dirichlet
h = L/partition[1]
γ = 10.0/h
a₁(u,v) = ∫(∇(v)⋅∇(u))dΩ +
          ∫( γ*(v*u) - (nl⋅∇(u))*v - (nl⋅∇(v))*u)dΓl +
          ∫( γ*(v*u) - (nr⋅∇(u))*v - (nr⋅∇(v))*u)dΓr +
          ∫( γ*(v*u) - (nb⋅∇(u))*v - (nb⋅∇(v))*u)dΓb +
          ∫( γ*(v*u) - (nt⋅∇(u))*v - (nt⋅∇(v))*u)dΓt
b₁(v) = ∫(f*v)dΩ +
        ∫( γ*(v*u₁) - (nl⋅∇(v))*u₁)dΓl +
        ∫( γ*(v*u₁) - (nr⋅∇(v))*u₁)dΓr +
        ∫( γ*(v*u₀) - (nb⋅∇(v))*u₀)dΓb +
        ∫( γ*(v*u₀) - (nt⋅∇(v))*u₀)dΓt
op₁ = AffineFEOperator(a₁,b₁,Uₕ,Vₕ)

# solve
uₕ₁ = solve(op₁)

# Output solution to vtk
writevtk(Ω,"solution",cellfields=["u_strong"=>uₕ₀, "u_weak"=>uₕ₁])

# Print error
println(√(∑(∫((uₕ₀-uₕ₁)*(uₕ₀-uₕ₁))dΩ)))


end
