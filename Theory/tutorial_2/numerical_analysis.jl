module Numerical_analysis
using Gridap

# Analytical solution
L = 3.0
EA = 1.0e3
uₑ(x) = x[1]^3*(L-x[1])^2
f(x) = -EA*Δ(uₑ)(x)

# Discrete model
model = CartesianDiscreteModel((0.0,L), (2,))
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"left",[1])
add_tag_from_tags!(labels,"right",[2])

# Triangulations
Ω = Triangulation(model)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(Ω,reffe,dirichlet_tags=["left","right"])
U = TrialFESpace(V,[uₑ,uₑ])

# Integration measures
dΩ = Measure(Ω,2)

# Weak form
a(u,v) = ∫(EA*(∇(v)⋅∇(u)))dΩ
l(v) = ∫(f*v)dΩ
op = AffineFEOperator(a,l,U,V)

# Solution
uₕ = solve(op)

# Consistency check
free_vals = rand(num_free_dofs(V))
dirichlet_vals = zeros(num_dirichlet_dofs(V))
vₕ = FEFunction(V,free_vals,dirichlet_vals)
println(abs(∑(a(uₑ,vₕ)-l(vₕ))))

# Galerkin orthogonality
println(abs(∑(a(uₕ-uₑ,vₕ))))

# L2 norm of the error
eₕ = uₕ-uₑ
l₂norm = √(∑(∫(eₕ*eₕ)dΩ))
println("L2 norm of the error: $(l₂norm)")

end
