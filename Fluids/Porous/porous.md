@def title = "Flow through porous zone"
@def mintoclevel=1
@def maxtoclevel=2 
@def hascode = true

<script>
  var coll = document.getElementsByClassName("collapsible");
  var i;

  for (i = 0; i < coll.length; i++) {
    coll[i].addEventListener("click", function() {
      this.classList.toggle("active");
      var content = this.nextElementSibling;
      if (content.style.display === "block") {
        content.style.display = "none";
      } else {
        content.style.display = "block";
      }
    });
  }
</script>

test

*Authors:*
*[Joël Ruesen](https://github.com/joeljulian) and [Oriol Colomés Gené](https://github.com/oriolcg)*

*Published:* January 2022

*Gridap version:* [Gridap@0.16.5](https://github.com/gridap/Gridap.jl/tree/release-0.14)

This tutorial shows how wave-progression through a porous medium is modelled. The model uses viscous incompressible Navier Stokes in combination with Darcy-Forchheimer resistance terms in the momentum balance, implemented using the Gridap library.

<!-- ~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/image.png" alt="Domain model" width="75%" /><center><i>3D model</i></center>~~~ -->

\toc


# Problem statement
# Mathematical formulation
Look at this shit
$$
% %\begin{equation}
\begin{split}
    \rho \left ( \frac{\partial \mathbf{u}}{\partial t} +\mathbf{u} \cdot \nabla \mathbf{u} \right ) -\mu \Delta \mathbf{u}+\nabla p=\mathbf{f} \quad \text{ on } \Omega \\
\nabla \cdot \mathbf{u}=0 \quad \text{ on } \Omega
\end{split}
% %\end{equation}
$$

## Governing equations
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
$$
% %\begin{equation}\label{eq:particle_velocity}
    \begin{split}
        u_{\text{hor}} (x,y,t) &= a_i \omega_i \frac{\cosh{(k_i y)}}{\sinh{(k_i d_{\text{water}})}} \sin(k_i x - \omega_i t - \theta_i)\\ 
        u_{\text{ver}} (x,y,t) &= a_i \omega_i \frac{\sinh{(k_i y)}}{\sinh{(k_i d_{\text{water}})}} \cos(k_i x - \omega_i t - \theta_i)
    \end{split}
% %\end{equation}
$$
### Free surface

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
$$
%\begin{equation}\label{eq:FS_BC_dyn}
\begin{split}
    \tau_{\text{fluid}} &= \tau_{\text{air}} \\
    -(\sigma_{\text{fluid}}) n_{\text{fa}} &= \sigma_{\text{air}} n_{\text{fa}} \\
    -(2 \mu \varepsilon(\mathbf{u}) - p I) n_{\text{fa}} &= (p_{\text{atm}} + \rho_{\text{fluid}} g \eta) n_{\text{fa}} \\
\end{split}
%\end{equation}
$$

### Sponge zone

# Numerical model
~~~<img style="display: block;max-width: 150%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/structure_code.PNG" alt="Domain model" width="100%" /><center><i>Structure of numerical components</i></center>~~~
## Input

## Weak form

### Weak Residual
```julia
#RESIDUAL
a((u,ut,p),(v,q))       = ∫( v⋅ut + v⋅((∇(u)')⋅u) + 2*ν*(ε(v)⊙ε(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ         # Original transient NS
b((u,p),(v,q))          = ∫( v⊙(αu)+ v⊙(β((u⋅u).^(1/2))*u) )dΩveg                                 # Darcy and forchheimer term
                          ∫( v⊙(C_damp*u))dΩdampout                                               # Sponge layer
c((u,p,η),(v,q,w))      = ∫( ( (-p_atm - g*η)*n_fa)⋅v )dΓ_fa                                       # traction equilibrium air/downwards
d((u,η,ηt),(v,w))       = ∫( ( (u⋅n_z) - ηt)*w )dΓ_fa                                              # traction equilibrium
l(v)                    = ∫( v⋅f_g )dΩ                                                             # Gravity throughout domain
e(v,u,p)                = ∫( v⋅(hN) )dΓ_out                                                        # Outflow boundary condition
tbt_oss((u,s),(v,ϵ))    = ∫( τₘ*((∇(u)')⋅u - s)⋅( (∇(v)')⋅u - ϵ ) )dΩ                               # Stabilization projection term

res(t,((u,s,p,η),(ut,st,pt,ηt)),(v,ϵ,q,w)) = a((u,ut,p),(v,q)) + b((u,p),(v,q)) - l(v) + e(v,u,p) - c((u,p,η),(v,q,w)) + d((u,η,ηt),(v,w)) + tbt_oss((u,s),(v,ϵ))
```
### Jacobian
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
## Domain
~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/domain.PNG" alt="Domain model" width="100%" /><center><i>Numerical wave tank</i></center>~~~
## Boundaries
## FE Spaces
### Triangulations
### Quadratures


# Solver
```julia
op = TransientFEOperator(res,jac,jac_t,X,Y)
xh0 = interpolate_everywhere([VectorValue(0.0,0.0),VectorValue(0.0,0.0),0.0,0.0],X(0.0))

ls = LinearSolvers.GmresSolver(verbose=false,preconditioner=ilu;τ=1e-6)
nl_solver = NLSolver(ls, show_trace = true,method = :newton,linesearch = BackTracking())
ode_scheme_1 = ThetaMethod(nl_solver, dt1, theta)
solver_1 = TransientFESolver(ode_scheme_1)

xh_1 = solve(solver_1, op, xh0, t0, T)
```

# Result analysis

The solution field can be saved in arrays of the variables over time:

```julia
# Setting up empty variables to push the values to
# Global variables to enable use in- and outside of for-loops
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
end
end # loop pushing data xh2
```

Furthermore, the simulation results can be saved in a JLD2 file to have all simulation results and their input parameters readily available together.

```julia
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
```

## Animations
For visualisation using ParaView, the package WriteVTK can be used

```julia
using WriteVTK

pvd = paraview_collection(datadir("sims", stage)*"\\$runname collection_ups", append=false)
pvd_Γ = paraview_collection(datadir("sims", stage)*"\\$runname collection_eta", append=false)

for (xh_tn, tn) in xh_1
        uh, sh, ph, ηh = xh_tn
        t_doc = round(tn; digits=3)
        pvd[tn] = createvtk(Ω, datadir("sims",stage)*"\\$runname FS$t_doc.vtu", cellfields = ["uh" => uh, "ph" => ph, "sh" => sh])
        pvd_Γ[tn] = createvtk(Γ_fa,datadir("sims",stage)*"\\$runname _FS_surf$t_doc.vtu",cellfields = ["etah" => ηh])
        xh0 = xh_tn
    end
```
## Plotting
~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/Fluids\Porous/img/Results_Hdomain.png" alt="Domain model" width="75%" /><center><i>Structure of framework components</i></center>~~~

# Workflow tools

To structure the work around setting up a model, running simulations, analysing data and collaborating on code through Github, the following tools are recommended.
## Research project structuring
The package 
```julia
using DrWatson
@quickactivate "PorousMediumFlow"

include(projectdir("src","Sim_PM.jl"))
```

## Simulating in batches

```julia
allparams = Dict(
    "stage"=>["numerical_methods"],
    "dt1"=> [ 0.005],
    "t0" => [ 0.0],
    "T"  => [ 0.01],###
    "dt2" => [0.01],
    "t0_2"  => [0.01],
    "T_2"  => [5.0], #30
    "a_in" => [0.015],  #m 0.05
    "T_in" => [1.0],
    "a_cons" => [ 0.15],# 0.2, 0.15, 0.1, 0.05],
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

function sim_prep(d::Dict)
    @unpack stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in a_cons, b_cons, C1, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength = d
    etahs = sim_run(stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in, a_cons, b_cons, C1, C2, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength)
    return true
end

for (i, d) in enumerate(dicts)
    f1 = sim_prep(d)
end
```

The julia file that runs the simulations can then make used of the created dictionaries of input parameters:
```julia
function sim_run(stage, dt1, t0, T, dt2, t0_2, T_2, a_in, T_in, theta_in, a_cons, b_cons, C1, C2, cellm, length_tank, depth_tank, length_damp, vegtop, vegheight, vegfront, veglength, inlet_x_loc, beachpor)

    ######################
    # Actual simulations #
    ######################

end


```

## Timing outputs

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

@@colbox-blue
:warning: ~~~<b>Warning:</b>~~~      
                     
This is warning box text
@@
