@def title = "Flow through porous zone"
@def mintoclevel=1
@def maxtoclevel=2 
@def hascode = true

*Authors:*
*[Joël Ruesen](https://github.com/joeljulian) and [Oriol Colomés Gené](https://github.com/oriolcg)*

*Published:* January 2022

*Gridap version:* [Gridap@0.16.5](https://github.com/gridap/Gridap.jl/tree/release-0.14)

This tutorial shows how wave-progression through a porous medium is modelled. The model uses viscous incompressible Navier Stokes in combination with Darcy-Forchheimer resistance terms in the momentum balance, implemented using the Gridap library.

<!-- ~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/image.png" alt="Domain model" width="75%" /><center><i>3D model</i></center>~~~ -->

\toc

# Problem statement
Large-scale offshore kelp farms have the beneficial side-effect that they can be used for coastal protection, increasing workability, and potentially reduce design requirements regarding fatigue life. All of this is possible through the wave damping potential of densely seeded grow lines or nets in seaweed cultivation areas.

Wave height reduction over a vegetated area is due to the dissipation of wave energy. Different methodologies exist to capture this problem in a model. Some models are able to capture even individual plant dynamics, they are computationally way too expensive for modelling on a farm scale. Other methods rely on calibrating a bulk drag coefficient, which captures the effects for a limited range of conditions and vegetation characteristics. 

The following model generalizes the complex geometry on a small scale by approaching the farm as a porous medium. Hereby, the effect of the kelp on incident waves can be modelled through a porous medium, described by the Darcy-Forchheimer resistance terms in the momentum equation.
# Mathematical formulation
Navier Stokes for visous incompressible flow, derived from the Cauchy momentum equation. Under the assumption of conservtion of mass and the resulting mass continuity equation this eventually reduces to:
$$
% %\begin{equation}
\begin{split}
    \rho \left ( \frac{\partial \mathbf{u}}{\partial t} +\mathbf{u} \cdot \nabla \mathbf{u} \right ) -\mu \Delta \mathbf{u}+\nabla p=\mathbf{f} \quad \text{ on } \Omega \\
\nabla \cdot \mathbf{u}=0 \quad \text{ on } \Omega
\end{split}
% %\end{equation}
$$

## Governing equations
The original Navier Stokes equations are elaborated with Darcy Forchheimer terms that are linear and quadratic with velocity:
$$
% %\begin{equation}
\begin{split}
     \underbrace{\frac{\partial \mathbf{u}}{\partial t}}_{\text{Velocity change}} + \underbrace{\mathbf{u} \cdot \nabla \mathbf{u}}_{\text{Convection}}  - \underbrace{\nu \Delta \mathbf{u}}_{\text{Diffusion}}+\underbrace{\nabla p^*}_{\text{Pressure gradient}}=  \mathbf{f}/\rho \quad \text{ on } \Omega \\
     \text{where } \quad \mathbf{f}/\rho = \underbrace{ \alpha \mathbf{u} + \beta \mathbf{u}^2 }_{\text{Darcy-Forchheimer resistance}}+ \underbrace{\mathbf{g}}_{\text{Gravity}}  \quad  \quad  \quad  \quad \\
     \nabla \cdot \mathbf{u}=0 \quad \text{ on } \Omega \\
\end{split}
% %\end{equation}
$$

## Boundary conditions
### Wave generation
Based upon linear wave theory, wave generation is provided through the following velocity components:
$$
% %\begin{equation}\label{eq:particle_velocity}
    \begin{split}
        u_{\text{hor}} (x,y,t) &= a_i \omega_i \frac{\cosh{(k_i y)}}{\sinh{(k_i d_{\text{water}})}} \sin(k_i x - \omega_i t - \theta_i)\\ 
        u_{\text{ver}} (x,y,t) &= a_i \omega_i \frac{\sinh{(k_i y)}}{\sinh{(k_i d_{\text{water}})}} \cos(k_i x - \omega_i t - \theta_i)
    \end{split}
% %\end{equation}
$$
### Free surface
The model remains single-phase and the free surface is defined through a linearized transpiration boundary condition, much like the one for Airy wave theory:
~~~<img style="display: block;max-width: 40%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/BC_fs.PNG" alt="Domain model" width="50%" /><center><i>Linearized transpiration boundary condition</i></center>~~~

$$
%\begin{equation}
    \begin{split}
    \sigma_\text{fluid} &= \tau_\text{s,v} - p I\\
    &= 2 \mu \varepsilon(\mathbf{u}) - p I
    \end{split}
%\end{equation}
$$
And for air
$$
%\begin{equation}
    \begin{split}
    \sigma_\text{air} &= p_{\text{atm}} + \rho_{\text{fluid}} g \eta
    \end{split}
%\end{equation}
$$
Plugging the stresses above into the traction equilibrium, the formulation is found that can be implemented in the weak form in the model.
$$
%\begin{equation}\label{eq:FS_BC_dyn}
\begin{split}
    \tau_{\text{fluid}} &= \tau_{\text{air}} \\
    -(\sigma_{\text{fluid}}) n_{\text{fa}} &= \sigma_{\text{air}} n_{\text{fa}} \\
    -(2 \mu \varepsilon(\mathbf{u}) - p I) n_{\text{fa}} &= (p_{\text{atm}} + \rho_{\text{fluid}} g \eta) n_{\text{fa}} \\
\end{split}
%\end{equation}
$$


# Numerical model
The mathematical equations posed above are captured in a numerical scheme to be able to simulate the resulting system.
The structure of the numerical method, with all components needed to arrive at a solution is presented below.

~~~<img style="display: block;max-width: 150%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/structure_code.PNG" alt="Domain model" width="100%" /><center><i>Structure of numerical components</i></center>~~~

Before setting up and code, all packages are called to include their functions in the following code.

```julia
using DrWatson
@quickactivate "PorousMediumFlow"
# Load packages
using Gridap
using Gridap.Geometry
using GridapODEs
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using LineSearches: BackTracking
using WriteVTK
using LinearAlgebra
using DelimitedFiles, FileIO, Dates
using JLD2, DataFrames, TimerOutputs
using IncompleteLU, IterativeSolvers
include("LinearSolvers.jl")
```
## Input variables
The constants that remain the same for each simulation and thus do not need to be user-defined are defined first.
Additionally, the variables that are direct functions of the input variables can be stated already.
```julia
    ## Establish constants:
    theta = 0.5      
    g = 9.81
    n_z = VectorValue(0.0,1.0)
    ρ = 1025.0
    p_atm = (1.013*10^5)/ρ
    ν =  0.001/ρ
    f_g = VectorValue(0.0, -g)

    λ_in = g*(T_in^2)/(2*π)
    k_in = 2*π/λ_in
    ω_in  = sqrt(g*k_in) 

    ## Geometry and meshing
    dampinback = length_damp
    dampoutfront = length_tank-length_damp
    dampoutback = length_tank
    n = cellm*length_tank
    m = (cellm)*depth_tank    
    h_cell = 1/cellm
    dampinfront = 0.0
    vegbottom = vegtop-vegheight
    vegback = vegfront+veglength

    ## Alternatively, dimensions can be parameterized
    # vegfront = round(λ_in, digits = 1, base = 5)
    # length_tank = round(10*λ_in, digits = 1, base =5)
    # veglength = round(4*λ_in, digits = 1, base =5)
    # length_damp = round(4*λ_in, digits = 1, base =5)
```
## Domain
The equations stated above are impemented onto the domain. The model essentially makes use of a stationary, evenly-spaced mesh, but _Gridap_ allows for more complex definitions as well.
```julia
    ## Domain and model set-up
    domain = (0,length_tank,0,depth_tank)
    partition = (n,m)
    model = CartesianDiscreteModel(domain,partition)
    writevtk(model,model_dir*"\\$runname _domain")

    k = 2
    Ω = Triangulation(model)
    degree = 2*k
    dΩ = Measure(Ω,degree)

    labels = get_face_labeling(model)
```
### Boundaries
The boundaries need to be defined (and tagged) in order to implement the boundary conditions on them later on.
 
```julia
    add_tag_from_tags!(labels,"diri_in",[3,7,])
    add_tag_from_tags!(labels,"freesurf_in",[3,])
    add_tag_from_tags!(labels,"free_surf",[3,4,6,])
    add_tag_from_tags!(labels,"neum_out",[4,8,])
    add_tag_from_tags!(labels,"diri_wall",[1,2,5])

    ## Prepare free surface discrete model
    bgface_to_mask = get_face_mask(labels,[6],1)
    Γface_to_bgface = findall(bgface_to_mask)
    model_Γ = BoundaryDiscreteModel(Polytope{1},model,Γface_to_bgface)
    writevtk(model_Γ,model_dir*"\\$(runname)_boundary")

    Γ_fa = Triangulation(model_Γ)
    dΓ_fa = Measure(Γ_fa,degree)
    n_fa = get_normal_vector(Γ_fa)

    Γ_out = BoundaryTriangulation(model,tags="neum_out")
    dΓ_out = Measure(Γ_out,degree)
    n_Γout = get_normal_vector(Γ_out)
```
### Zones within domain
For the damping zone and actual porous medium, where the porous parameters are activated, the coordinates must be specified:

```julia
    function is_inveg(coords)
        midx = (coords[1][1] + coords[2][1])/2
        midy = (coords[1][2] + coords[3][2])/2
        midx >= vegfront && midx <= vegback && midy <= vegtop && midy >= vegbottom
    end

    function is_indampout(coords)
        midx = (coords[1][1] + coords[2][1])/2
        midy = (coords[1][2] + coords[3][2])/2
        midx >= dampoutfront && midx <= dampoutback
    end

    oldcell_to_coods = get_cell_coordinates(Ω)
    oldcell_to_is_inveg = [lazy_map(is_inveg, oldcell_to_coods)[i] for i in 1:length(lazy_map(is_inveg, oldcell_to_coods))]
    oldcell_to_is_indampout = [lazy_map(is_indampout, oldcell_to_coods)[i] for i in 1:length(lazy_map(is_indampout, oldcell_to_coods))]
    incell_to_cellveg = findall(oldcell_to_is_inveg)
    incell_to_celldampout = findall(oldcell_to_is_indampout)
    outcell_to_cell = findall(iszero,(oldcell_to_is_inveg+ oldcell_to_is_indampin+ oldcell_to_is_indampout))

    model_veg = DiscreteModel(model,incell_to_cellveg)
    model_dampout = DiscreteModel(model,incell_to_celldampout)

    Ωveg = Triangulation(model_veg)
    dΩveg = Measure(Ωveg,degree)
    Ωdampout = Triangulation(model_dampout)
    dΩdampout = Measure(Ωdampout,degree)

    # Vegetation parameters (currently user-defined in input)
    α_veg = a_cons
    β_veg = b_cons

    # Sponge damping zone
    α_dampout(x::VectorValue) = C1*(1-((depth_tank-x[2])-depth_tank)/(depth_tank))*((x[1]-(length_tank-length_damp))/(length_tank-length_damp))^2
```
### Inflow conditions
As stated, wave generation is implemented by forcing particle velocity at the leftmost boundary. This is done using the orbital velocities as defined using Linear wave theory.
```julia
    # #### Boundary conditions and porous zone
    # Velocity profile inlet
    u_hor(x,t) = a_in*ω_in*((cosh(k_in*x[2]))/(sinh(k_in*depth_tank)))*sin(x[1]*k_in - ω_in*t - theta_in)
    u_ver(x,t) = a_in*ω_in*((sinh(k_in*x[2]))/(sinh(k_in*depth_tank)))*cos(x[1]*k_in - ω_in*t - theta_in)
   
    # Boundary conditions at inlet and wall
    u_in(x, t::Real) = VectorValue(u_hor(x,t), u_ver(x,t))
    u_wall(x, t::Real) = VectorValue(0.0, 0.0)
    u₀(x, t::Real) = VectorValue(0.0,0.0)
    #u₀(x, t::Real) = VectorValue(a_in*ω_in*((cosh(k_in*x[2]))/(sinh(k_in*depth_tank))), 0.0)        #Initial conditions can be chosen

    u_in(t::Real) = x -> u_in(x, t)
    u_wall(t::Real) = x -> u_wall(x, t)
    u₀(t::Real) = x -> u₀(x, t::Real)
```
## FE Spaces
The spaces are defined on which the functions will need to be found that provide the eventual solution fields.
This is in line with the [Gridap tutorials](https://gridap.github.io/Tutorials/dev/). 
```julia
    # ## Numerical Scheme
    # Velocity FE space
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},k)
    V₀ = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["diri_wall", "diri_in"])
    U = TransientTrialFESpace(V₀, [u_wall, u_in])

    # Stabilization variable Space
    reffe_s = ReferenceFE(lagrangian,VectorValue{2,Float64},k)
    S₀ = TestFESpace(model, reffe_s, conformity=:H1)
    S = TrialFESpace(S₀)

    # Pressure FE space
    reffe_p = ReferenceFE(lagrangian,Float64,1)
    Q = FESpace(model, reffe_p, conformity=:C0)
    P = TrialFESpace(Q)

    # Free surface space
    reffe_η = ReferenceFE(lagrangian,Float64,k)
    Η₀ = TestFESpace(model_Γ, reffe_η, conformity=:H1) 
    Η = TrialFESpace(Η₀)

    # Combining in multi-field
    Y = MultiFieldFESpace([V₀,S₀,Q,Η₀])
    X = TransientMultiFieldFESpace([U,S,P,Η])
```
## Weak form

Before going on to implementing the physics, some auxiliary functions are defined:
```julia
    @timeit timesave "weak form setup" begin
    nl_floor(u) = max((u⋅u).^(1/2), 1e-10)
    c₁    = 12.0*ρ           # algorithmic constant (John,Kindl2007)
    c₂    = 2.0
    ca    = c₁*ν/(h_cell^2)  # aux. constants
    cb    = c₂/h_cell
    τₘ(u)       = 1/( ca + cb*((u⋅u).^(1/2)) )                                           # Stabilization parameter (Colomes2016)
    τlin(u,du)  = -((cb)/( (ca + cb*((u⋅u)^(1/2)) )^2)) * ((1/(nl_floor(u)))*(u⋅du))      # Linearized stabilization parameter # u⋅du/nl_floor(u)
    hN(x)       = VectorValue(p_atm + g*(depth_tank-x[2]), 0.0)                          # hydrostatic pressure
```
### Weak Residual
The governing equations are captured through the weak form:
```julia
    #RESIDUAL
    a((u,ut,p),(v,q))       = ∫( v⋅ut + v⋅((∇(u)')⋅u) + 2*ν*(ε(v)⊙ε(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ         # Original transient NS
    b((u,p),(v,q))          = ∫( v⊙(αu)+ v⊙(β((u⋅u).^(1/2))*u) )dΩveg                                 # Darcy and forchheimer term
                              ∫( v⊙(C_damp*u))dΩdampout                                               # Sponge layer
    c((u,p,η),(v,q,w))      = ∫( ( (-p_atm - g*η)*n_fa)⋅v )dΓ_fa                                       # traction equilibrium air/downwards
    d((u,η,ηt),(v,w))       = ∫( ( (u⋅n_z) - ηt)*w )dΓ_fa                                              # traction equilibrium
    l(v)                    = ∫( v⋅f_g )dΩ                                                             # Gravity throughout domain
    e(v,u,p)                = ∫( v⋅(hN) )dΓ_out                                                        # 'Outflow' pressure profile
    tbt_oss((u,s),(v,ϵ))    = ∫( τₘ*((∇(u)')⋅u - s)⋅( (∇(v)')⋅u - ϵ ) )dΩ                               # Stabilization projection term

    res(t,((u,s,p,η),(ut,st,pt,ηt)),(v,ϵ,q,w)) = a((u,ut,p),(v,q)) + b((u,p),(v,q)) - l(v) + e(v,u,p) - c((u,p,η),(v,q,w)) + d((u,η,ηt),(v,w)) + tbt_oss((u,s),(v,ϵ))
```
### Jacobian
Since Backtracking line search will be used in the solver, the derivative of the weak form is also needed; the jacobian:
``` julia
    #JACOBIAN
    jac_a(u,(du,dp),(v,q))      = ∫( v⋅((∇(du)')⋅u) + v⋅((∇(u)')⋅du) + 2*ν*(ε(v)⊙ε(du)) - (∇⋅v)*dp + q*(∇⋅du))dΩ
    jac_b(u,du,v)               = ∫( v⊙(αdu) + (β*v⋅((1/(nl_floor∘(u)))(u⊗u)⋅du + (u⋅u).^(1/2)*du)) )dΩveg +
                                  ∫( v⊙(C_damp*du)  )dΩdampout
    jac_c(dη,v)                 = ∫( ((-g*dη)*n_fa)⋅v )dΓ_fa
    jac_d(ηt,du,w)              = ∫( ( (du⋅n_z) - ηt)*w )dΓ_fa
    jac_tbt_oss(u,s,du,ds,v,ϵ)  = ∫( τₘ*( ((∇(v)')⋅du)⋅((∇(u)')⋅u - s) + (((∇(v)')⋅u-ϵ)⋅((∇(u)')⋅du +(∇(du)')⋅u -ds)))
                                + (τlin∘(u,du))*(((∇(u)')⋅u - s)⋅( (∇(v)')⋅u - ϵ ) ))dΩ

    jac(t,((u,s,p,η),(ut,st,pt,ηt)),(du,ds,dp,dη),(v,ϵ,q,w))        = jac_a(u,(du,dp),(v,q))  + jac_b(u,du,v) + jac_c(dη,v) + jac_d(ηt,du,w) + jac_tbt_oss(u,s,du,ds,v,ϵ) #+ ∫( v⋅(dp*n_Γout-2*ν*ε(du)⋅n_Γout) )dΓ_out
    jac_t(t,((u,s,p,η),(ut,st,pt,ηt)),(dut,dst,dpt,dηt),(v,ϵ,q,w))  = ∫( v⋅dut )dΩ - ∫( dηt*w )dΓ_fa
```
The resulting system provides the operator:
```julia
    op = TransientFEOperator(res,jac,jac_t,X,Y)
```
### Initial conditions
```julia
    # ###### Stokes operator
    res_stokes((u,p),(v,q)) = ∫( 2*ν*(ε(v)⊙ε(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ
    jac_stokes((u,p),(du,dp),(v,q)) = res_stokes((du,dp),(v,q))
    op_stokes = FEOperator(res_stokes,jac_stokes,X(0.0),Y(0.0))
    #global xh0 = solve(op_stokes)
    xh0 = interpolate_everywhere([VectorValue(0.0,0.0),VectorValue(0.0,0.0),0.0,0.0],X(0.0))
    #xh0 = interpolate_everywhere([u₀(0.0),VectorValue(0.0,0.0),0.0,0.0],X(0.0))
    #writevtk(Ω,datadir("sims", "$runname")*"\\$runname _initial_up",cellfields=["u"=>xh0[1],"p"=>xh0[2]])
```
## Solver

```julia
    # ## Solver
    ls = LinearSolvers.GmresSolver(verbose=false,preconditioner=ilu;τ=1e-6)
    nl_solver = NLSolver(ls, show_trace = true,method = :newton,linesearch = BackTracking())
    
    ## Without preconditioning
    # nl_solver = NLSolver(show_trace = true,method = :newton,linesearch = BackTracking())

    ode_scheme_1 = ThetaMethod(nl_solver, dt1, theta)
    solver_1 = TransientFESolver(ode_scheme_1)
    xh_1 = solve(solver_1, op, xh0, t0, T)
    
    pvd = paraview_collection(datadir("sims", stage)*"\\$runname collection_ups", append=false)
    pvd_Γ = paraview_collection(datadir("sims", stage)*"\\$runname collection_eta", append=false)

    for (xh_tn, tn) in xh_1
        uh, sh, ph, ηh = xh_tn
        t_doc = round(tn; digits=3)
        pvd[tn] = createvtk(Ω, datadir("sims",stage)*"\\$runname FS$t_doc.vtu", cellfields = ["uh" => uh, "ph" => ph, "sh" => sh]) #
        pvd_Γ[tn] = createvtk(Γ_fa,datadir("sims",stage)*"\\$runname _FS_surf$t_doc.vtu",cellfields = ["etah" => ηh]) #
        # global xh0
        # xh0 = xh_tn
    end
```
The code above will produce a solution that's immediately written to VTK files.
By writing the solution to variable `xh0`, this can now be loaded as the intial conditions for a new solution.

Essentially this provides the option to first run a small timespan with large timesteps to provide a rough solution.
The following run can then performed at higher temporal resolution, while it can skip the (uninteresting) initial start-up behaviour.

```julia
    pvd = paraview_collection(datadir("sims",stage)*"\\collection_ups_$(runname)", append=false)
    pvd_Γ = paraview_collection(datadir("sims",stage)*"\\collection_eta_$(runname)", append=false)

    # Same solver method, new inputs for Phase 2
    ode_scheme_2 = ThetaMethod(nl_solver, dt2, theta)
    solver_2 = TransientFESolver(ode_scheme_2)
    xh_2 = solve(solver_2, op, xh0, t0_2, T_2)

    global ηns = []
    global uhs1 = []
    global uhs2 = []
    global uhs3 = []
    global uhs4 = []
    global phs = []
    for (i, (xh, t)) in enumerate(xh_2)
        uh, sh, ph, ηh = xh
    
        global cell_values_ηn = get_cell_dof_values(ηh)
        surface = []
        len_FS = length(cell_values_ηn)
            for j in 1:len_FS
                push!(surface, cell_values_ηn[j][3])
            end
        push!(ηns, surface')
    
        global cell_values_uh = get_cell_dof_values(uh)
        surface_uh1 = []
        surface_uh2 = []
        surface_uh3 = []
        surface_uh4 =[]
        len_UH = length(cell_values_uh)
            for j in 1:len_UH
                push!(surface_uh1, cell_values_uh[j][1])
                push!(surface_uh2, cell_values_uh[j][2])
                push!(surface_uh3, cell_values_uh[j][3])
                push!(surface_uh4, cell_values_uh[j][4])
            end
        push!(uhs1, surface_uh1')
        push!(uhs2, surface_uh2')
        push!(uhs3, surface_uh3')
        push!(uhs4, surface_uh4')
    
        global cell_values_ph = get_cell_dof_values(ph)
        surface_ph = []
        len_pH = length(cell_values_ph)
            for j in 1:len_pH
                push!(surface_ph, cell_values_ph[j][1])
            end
        push!(phs, surface_ph')
    
        t_doc = round(t; digits=4)
        pvd[t] = createvtk(Ω, datadir("sims",stage)*"\\$runname _FS$t_doc.vtu", cellfields = ["uh" => uh, "ph" => ph, "sh" => sh])
        pvd_Γ[t] = createvtk(Γ_fa,datadir("sims",stage)*"\\$runname _FS_surf$t_doc.vtu",cellfields = ["etah" => ηh])
    end
    
    vtk_save(pvd)
    vtk_save(pvd_Γ)
```
We now have the desired solution as variable `xh`, as timeseries in the variables `ηns, uhs1, uhs2, uhs3, uhs4` and `phs`, and as the saved VTK files.
To be able to process the solution at a later moment as wel, JLD2 files are created. Native julia datafiles that allow for storing the generated solution. To this end, the arrays are written into a library that is then saved in the `.jld2` format. 

```julia
    @timeit timesave "Data build & save" begin
    
    dat_total = Dict()
    dat_total["eta"] = ηns
    dat_total["uh1"] = uhs1
    dat_total["uh2"] = uhs2
    dat_total["uh3"] = uhs3
    dat_total["uh4"] = uhs4
    dat_total["ph"] = phs
    dat_params = Dict(
        "stage" => [stage],
        "runname" => [runname],
        "rundate" => [start_time],
        "Lwave"=> [λ_in],
        "dt1"=> [ dt1],
        "t0" => [ t0],
        "T"  => [T],
        "dt2" => [dt2],
        "t0_2"  => [t0_2],
        "T_2"  => [T_2],
        "a_in" => [ a_in],
        "T_in" => [ T_in],
        "theta_in" => [theta_in],
        "a_cons" => [a_cons],
        "b_cons" => [b_cons],
        "C1" => [ C1],
        "cellm"  => [cellm],
        "length_tank"  => [length_tank],
        "depth_tank"  => [depth_tank],
        "length_damp" => [length_damp], 
        "vegtop" => [vegtop], 
        "vegbottom" => [vegbottom],
        "vegheight" => [vegheight], 
        "vegfront" => [vegfront], 
        "vegback" => [vegback], 
        "veglength"=>[veglength]
        )
    dat_total["params"] = dat_params

    save(datadir("data",stage)*"\\$(runname).jld2",dat_total)

    println("Starttime was: $(start_time), current time is ", Dates.format(Dates.now(),"HH_MM_SS"))
    return true
end
end
```

## Plotting
Julia has many options for graphically presenting results. Most are julia-versions of software you may have used before (think matplotlib).
[Plots](http://docs.juliaplots.org/latest/) is a visualization interface and toolset. It sits above other backends, like GR, PyPlot, PGFPlotsX, or Plotly.

For the numerical model, a series of functions was defined to plot various insights. An indication of the core functions is provided here.
```julia
using Plots; gr()
function plot_wave(loadedruns,region,starttime,saveorshow,plotaddition)
    global params1 = loadedruns[1]
    plotname = "Damping_$(plotaddition)"
    titlename = "Reproduction - Zhu et al. (2021), waveAT015T2"
    xaxisname = "Normalized distance into vegetation (x-Lₛ)/Lᵥ"
    yaxisname = "Normalized Wave height η/η₀"
    yaxisname = "Wave height η [m]"
    labelsname = ["Hs = 1.5 m ,T = 6 s, α = $(params1["params"]["a_cons"][1])"]#, β = $(params1["params"]["b_cons"][1])"]#,"T = $(params2["params"]["T_in"][1]), a = $(params2["params"]["a_cons"][1]), b = $(params2["params"]["b_cons"][1])","T = $(params3["params"]["T_in"][1]), a = $(params3["params"]["a_cons"][1]), b = $(params3["params"]["b_cons"][1])",
    labelsname = ["a_in = $(params1["params"]["a_in"][1]),T = $(params1["params"]["T_in"][1]), α = $(params1["params"]["a_cons"][1]), β = $(params1["params"]["b_cons"][1])"]#,"T = $(params2["params"]["T_in"][1]), a = $(params2["params"]["a_cons"][1]), b = $(params2["params"]["b_cons"][1])","T = $(params3["params"]["T_in"][1]), a = $(params3["params"]["a_cons"][1]), b = $(params3["params"]["b_cons"][1])",
    global max_at_coord1 = []
    
    global finalcell = []
    global cellm_list = []
    global cphase = (params1["params"]["Lwave"][1]/params1["params"]["T_in"][1])
```
```julia
    if starttime ==1000
        temp_time = 2*trunc(Int,params1["params"]["length_tank"][1]/cphase)
        if temp_time >= params1["params"]["T_2"][1]
            println("initial wave travelling time is $(temp_time)s, choosing last 2s of simulation")
            global start_time = 2*trunc(Int,params1["params"]["length_tank"][1]/cphase)-1
        else
        global start_time = 2*trunc(Int,params1["params"]["length_tank"][1]/cphase)
        println("No measurement starttime inserted: using c=$(cphase)m/s -> t=$(start_time)")
        end
    elseif starttime ==0
        global start_time = 1
    else
        global start_time = starttime
    end

    if length(params1) == 0
        println("no data loaded yet, commencing now")
        global params1 = load(datadir("data",foldername,runs[1]))
        # global xcoord_cell = trunc(Int,(params1["params"]["vegfront"][1]+params1["params"]["veglength"][1])*params1["params"]["cellm"][1])
    else
        println("data already loaded, skipping step")
    end
        println("Houston, we have data")
   
        # for i in 1:length(runs)
        # global run_dict = load(datadir("data",foldername,runs[i]))
        # icellm = run_dict["params"]["cellm"][1]
        # push!(cellm_list, icellm)
        # ifinalcell = run_dict["params"]["length_tank"][1]*icellm
        # push!(finalcell, ifinalcell)
```
Depending on what type of plot is desired, the entire domain, or just part of it can be shown.
```julia
    if region == "veg"
        startcell = trunc(Int,params1["params"]["vegfront"][1]*params1["params"]["cellm"][1])
        endcell = trunc(Int,startcell + params1["params"]["veglength"][1]*params1["params"]["cellm"][1])
        totalcells = endcell-startcell        
        EDRfirstcell = startcell
        EDRlastcell = totalcells
        x_axis = LinRange(params1["params"]["vegfront"][1], (params1["params"]["vegfront"][1]+ params1["params"]["veglength"][1]), totalcells)
        k=1
        xaxisname = "Normalized distance into vegetation (x-Lₛ)/Lᵥ"

    elseif region == "domain"
        startcell = 3 #was1
        endcell = trunc(Int,params1["params"]["length_tank"][1]*params1["params"]["cellm"][1])
        EDRfirstcell = trunc(Int, params1["params"]["vegfront"][1]*params1["params"]["cellm"][1])
        EDRlastcell = EDRfirstcell + trunc(Int,(params1["params"]["veglength"][1]*params1["params"]["cellm"][1]))
        totalcells = endcell-startcell
        x_axis = LinRange(0,params1["params"]["length_tank"][1],totalcells)
        # x_axis = LinRange(0,750,totalcells)

        k=4
        xaxisname = "Distance into domain [m]"

    else 
        printlnln("No region selected")
    end
```
Next up is looping through the timesteps to build the arrays that will be plotted eventually.
```julia
        for j in 1:totalcells
            global timeseries = []
            for ts in trunc(Int,(start_time)/params1["params"]["dt2"][1]):(trunc(Int,(params1["params"]["T_2"][1])/params1["params"]["dt2"][1])-1)
                temp_val = params1["eta"][ts][j+startcell]
                push!(timeseries, temp_val')
            end
            temp_max_at_coord = maximum(timeseries)
                push!(max_at_coord1, temp_max_at_coord)
        end
    # end

    norm_run1 = (max_at_coord1/max_at_coord1[k])
    
    #If desired, a linear fit can be added:
    linearizedfit, coef1, coef2 = poly_fit(x_axis[EDRfirstcell:EDRlastcell],max_at_coord1[EDRfirstcell:EDRlastcell])
    coefn1, coefn2 = [coef1, coef2] ./ max_at_coord1[k]
    linearizedfitnorm = linearizedfit/max_at_coord1[k]
    
    MP1, MP2 = x_axis[EDRlastcell], linearizedfitnorm[EDRlastcell-EDRfirstcell]
    perc1 = round((1-(linearizedfit[EDRlastcell-EDRfirstcell]/max_at_coord1[k])^2)*100, digits = 1)
    HTRend = round(linearizedfit[EDRlastcell-EDRfirstcell]/max_at_coord1[k], digits = 2)
    # lenslo = 0.95*minimum([norm_run1[end], norm_run2[end], norm_run3[end], norm_run5[end], norm_run5[end], norm_run6[end]])#, norm_run4[end]])
    # lenshi = 1.05*maximum([norm_run1[end], norm_run2[end], norm_run3[end], norm_run5[end], norm_run5[end], norm_run6[end]])#, norm_run4[end]])
    # lensleft = 7.5/8
    # lensright = 1
    # xlens = 0.4
    # ylens = 0.3
    # lenswidth = 0.2
    # lensheight = 0.35

    global labels = labelsname    
    global p_total1 = plot(x_axis, norm_run1, label = labels[1], dpi=200, legend = :bottomleft, ylim=(0,1.3))
    # plot!(x_axis[EDRfirstcell:EDRlastcell], linearizedfitnorm, ls=:dash, label="linear fit, lsq ($(trunc(coefn1,digits=3)) + $(trunc(coefn2, digits=3))x)")
    scatter!(xpoint_exp_wave, values_exp_wave, label = "Data Hu et al. (2014)")
    scatter!([MP1], [MP2], markershape=:star5, color=:grey, alpha=0.5, label = "HTR = $(HTRend), EDR = $(perc1) %")
    xlabel!(p_total1,xaxisname)
    ylabel!(p_total1, yaxisname)
    title!(p_total1, titlename)
    # lens!([lensleft,lensright], [lenslo, lenshi], inset = (1, bbox(xlens, ylens, lenswidth, lensheight)), xticks = false, framestyle=:box, ls=:dash, lc=:grey, subplot = 2)#,xticks=1:2:2, subplot=1)

    if saveorshow == "save"
        savefig(p_total1, plotsdir("$(plotname).png"))       
        println("figure saved as $( plotsdir("$(plotname).png"))")
    elseif saveorshow == "show"
        println("not saved")
    else
        println("please define show or save")
    end
    p_total1
end
```
### Fitting a polynomial
For the analysis of data and for further development of the framework, it is useful to be able to fit user-specified function types to the simulation data. This can be done using many packages, [LsqFit](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) is one convenient option.

```julia
using Polynomials
using LsqFit
function poly_fit(xdata,ydata)

    p0 = [0.1, 0.1]
    model(t,p)= p[1] .+ p[2]*t
    xlin = range(xdata[1], xdata[end], length =100)

    fitcurve = curve_fit(model, xdata, ydata, p0)
    residuals = fitcurve.resid

    global linydata = ydata .+ residuals
    param = fitcurve.param
    println("$(param)")

    return linydata, param[1], param[2]
end
```

~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/Results_Hdomain.png" alt="Domain model" width="75%" /><center><i>Structure of framework components</i></center>~~~

# Workflow tools

To structure the work around setting up a model, running simulations, analysing data and collaborating on code through Github, the following tools are recommended.
## Research project structuring
The package [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/dev/) automatically creates a folder-structure/repository that includes everything you need to develop code in a research project and allows for collaboration through GitHub easily. The set-up includes infrastructure dedicated to ignore some folders during commits to avoid overloading your repo. It lets anyone duplicate all package versions and dependencies for reproducing your work without running into errors. Furthermore, provide clear naming conventions and prevent double work by checking whether the same input parameters have already been used in an earlier run.

It is worth checking out and investing a little bit of time in beforehand, as it will make your life significantly easier.
```julia
using DrWatson
@quickactivate "PorousMediumFlow" # defines what environment (ie. versions etc) is activated
include(projectdir("src","Sim_PM.jl")) #projectdir() is an example of DrWatson's helpful tools
```

## Simulating in batches
With simulations that take more than a couple of seconds, automation is your friend.
DrWatson allows to setup a library of input parameters that are fed into the simulation script when the last run is done.

```julia
allparams = Dict(
    "stage"=>["numerical_methods"],
    "dt1"=> [ 0.005],   #timestep of first solving round
    "t0" => [ 0.0],     #starttime of first solving round
    "T"  => [ 0.01],    #end time of first solving round

    "dt2" => [0.01],
    "t0_2"  => [0.01],
    "T_2"  => [5.0], #30

    "a_in" => [0.015, 0.01],  #m 0.05
    "T_in" => [1.0],
    
    "a_cons" => [ 0.15 0.2, 0.1, 0.05],
    "b_cons" => [ 5.0],
    "C1" => [ 7.5],

    "cellm"  => [5.0],
    "length_tank"  => [15.0],
    "depth_tank"  => [1.0],
    "length_damp" => [7.0], 
    "vegtop" => [1.0], 
    "vegheight" => [0.1], 
    "vegfront" => [1.0], 
    "veglength" => [6.0]
)
dicts = dict_list(allparams)
```
We have now set up the dictionary that contains all variables we desire for input.
We run `sim_prep(d:Dict)` to build a series of _all_ combinations of these variables. 

@@colbox-blue
:warning: ~~~<b>Warning:</b>~~~      

For the large number of input variables of this model, the possible combinations grow quickly.
Make sure to carefully think about what runs (and combinations) are essential.
@@

```julia
function sim_prep(d::Dict)
    @unpack stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in a_cons, b_cons, C1, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength = d
    etahs = sim_run(stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in, a_cons, b_cons, C1, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength)
    return true
end

for (i, d) in enumerate(dicts)
    f1 = sim_prep(d)
end
```

The julia file that runs the simulations can then make used of the created dictionaries of input parameters:
```julia
function sim_run(stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in, a_cons, b_cons, C1, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength)
    ######################
    # Actual simulations #
    ######################
end
```

## Timing outputs
[TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) is a minimal, basic tool that helps to gain some insights into the resources your simulations demands.

```julia
using TimerOutputs

timesave = TimerOutput()
@timeit timesave "Total simulation time" begin

    ####################################
    # Setting up a batch of simulation #
    ####################################

    @timeit timesave "Simulation time 1 run" begin

        ######################
        # Actual simulations #
        ######################

    end
end
show(timesave)
```
The image below provides an example of what the results look like when a .jl uses multiple timers throughout:

~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/base_run.png" alt="TimerOutputs summary" width="100%" /><center><i>TimerOutputs summary</i></center>~~~

