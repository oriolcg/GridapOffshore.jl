@def title = "Very Large Floating Structure (VLFS)"
@def mintoclevel=1
@def maxtoclevel=2 
@def hascode = true

*Authors:*
*[Dorette Regout](https://github.com/DoretteR) and [Oriol Colomés Gené](https://github.com/oriolcg)*

*Published:* October 2021

*Gridap version:* [Gridap@0.14](https://github.com/gridap/Gridap.jl/tree/release-0.16)

This tutorial shows how to solve a Fluid Structure Interaction (FSI) problem using Gridap and provides the instructions to build a 2D model considering a multi-module VLFS, solved in the frequency domain.   

\toc

# Problem statement

As shown in the figure below, the multi-module VLFS is represented by four floating one-dimensional beams, interconnected by rotational springs, where the structure is located on top of a two-dimensional fluid domain. The dynamic response of the structure is described by the Euler-Bernoulli beam theory and it is assumed that the beams will only perform small displacements in vertical direction. As the structure is subject to small, single frequency, incident head waves (with a wavelength $\lambda$), the linear (Airy) wave theory is therefore assumed to be valid, provided that the beams are always in contact with the fluid.

The total length of the VLFS is denoted by $L$, where the beams have a uniform and homogeneous cross-section along the structure. The bending stiffness, the mass density, the height, and the length of each beam element are denoted by $EI$, $\rho_b$, $h_b$, and $βL$, respectively. The rotational stiffness for each connection is defined by $k_r = ξEI/L$, where $\xi$ is the rotational stiffness parameter. The 2D fluid domain, denoted as $\Omega$, has a constant depth, $d$, with a constant density, $\rho_w$. The boundaries of the domain are defined by the seabed, $\Gamma_b$, the free water surface, $\Gamma_{fs}$, the fluid-structure interface, $\Gamma_{str}$, and two vertical boundaries, $\Gamma_{L}$ and $\Gamma_{R}$, on the left and right side of the domain. The fluid is considered to be inviscid, incompressible, and has irrotational flow. As a result, the fluid is
expressed by the velocity potential, which is a scalar function that represents the velocity field of the water in terms of spatial derivatives of the scalar function $\phi$.  


~~~<img style="display: block;max-width: 100%;height: auto;margin: auto;float: none!important;" src="/FSI/VLFS/img/numerical_model.png" alt="VLFS model" width="75%" /><center><i>VLFS model</i></center>~~~


# Transformation to the frequency domain

Due to the excitation of the waves, the water pressure and the elevation at the water surface will change, causing vertical displacements of the floating beams. This interaction between the motion of the water and the displacement of the structure is mathematically described by partial differential equations (PDEs) that initially depend on space and time. Accordingly, initial conditions (with respect to time) and boundary/interface conditions (with respect to space) have to be determined to solve the problem. In another tutorial called ''[Very Flexible Floating Structure](/FSI/VFFS)'', it is demonstrated how to solve the FSI problem for a single floating 1D beam in the time-domain.

However, solving these type of FSI problems in the time domain may be computationally intensive and it is therefore more convenient to solve the problem in the frequency domain. In general, time dependent PDEs can be transformed to the corresponding frequency domain by applying the following [Fourier transform pair] (https://en.wikipedia.org/wiki/Fourier_transform):

$$
\begin{equation*}
    \tilde{G}(\omega) = \int_{-\infty}^{\infty} g(t)\,e^{-i\omega t}\,\,dt \qquad \text{and} \qquad  g(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{G}(\omega)\,e^{i\omega t}\,\,d\omega
\end{equation*}
$$


Due to the linearity of the system, this approach may be applied and each quantity in the time domain is rewritten in its respective term in the frequency domain. Hence, all equations for this model are expressed in terms of their space dependence in the frequency domain, where all quantities are complex-valued and are frequency dependent. Consequently, the solution to the problem will be solved regarding the steady-state solution.




# Mathematical formulation
To build the FSI model, first, the mathematical formulation to describe the interaction of the coupled system is discussed. This results a system of PDEs, where this set of equations will include the following components:
- equation of motion (EoM) of the fluid
- equation of motion (EoM) of the structure
- boundary/interface conditions (BCs/ICs) for the fluid domain
- boundary/interface conditions (BCs/ICs) for the structural domain


## EoM fluid
Based on the aforementioned fluid characteristics, the divergence of the fluid velocity is considered to be zero. This results in the velocity potential $\phi$ to satisfy the Laplace equation. Similarly, the fluid pressure can be expressed in terms of the velocity potential, where the equation for the pressure is derived from the momentum balance and the linearized Bernoulli equation (based on the linear wave theory):


$$
\left\{\begin{array}{l}
\nabla \cdot \vec{u}=0 \\
\nabla \phi=\vec{u}
\end{array} \Leftrightarrow \nabla \cdot(\nabla \phi)=\Delta \phi=0 \quad\right. \text { in } \quad \Omega
$$


$$
\begin{align*}
    &-i\omega\phi + \frac{p}{\rho_w} + gz = 0 \qquad \quad \qquad \text { in } \quad \Omega\\
    \rightarrow \quad &p \,=\, -\rho_w gz \, +\, iω ρ_w  ϕ  
\end{align*}
$$


## EoM VLFS
According to the Euler-Bernoulli beam theory, the equation of motion of the VLFS is described as follows:

$$
\begin{equation*}
    -\omega^2\rho_{b}h_b w + EI\frac{\partial^4w}{\partial x^4} = p\,\big|_{\Gamma_{str}} \quad \quad \text{on}\quad \Gamma_{str}
\end{equation*}
$$

Where $w$ is the vertical displacement of the VLFS in z-direction and $p$ is the pressure at the water surface, acting as a distributed load on the bottom of the VLFS, i.e. the interface of the structure and the fluid.


To solve the problem, the following boundary conditions for both the fluid and the VLFS structure are determined:

## Sea bed boundary
The impermeability condition holds at the seabed, which means that no water is allowed to flow through the bottom. Therefore, the flow velocity normal to the seabed is equal to zero. As the velocity field of the water can be expressed in terms of spatial derivatives of the scalar function $\phi$, this results in:

$$
\begin{equation*}
    \vec{n}\cdot\nabla\phi = 0 \qquad \text{on} \quad \Gamma_b
\end{equation*}
$$


## Free surface
At the free surface, the water should satisfy both a kinematic and a dynamic boundary condition. To this end, a function for the surface elevation of the water is introduced; where $η$ is the elevation in water surface, measured relative to the mean water level.

For the kinematic boundary condition holds that the water particles cannot 'leave' the water surface. Therefore, the flow velocity in normal direction to the surface should be equal to the velocity of the surface elevation. Regarding the dynamic condition, the water pressure at the surface should be equal to the atmospheric pressure, i.e. equal to zero:

$$
\begin{align*}
    \text{kinematic:}& \qquad \vec{n}\cdot\nabla\phi = -iωη  &&\text{on} \quad \Gamma_{fs}\\
    \text{dynamic:}& \qquad -i\omega\phi + g\eta = 0 &&\text{on} \quad \Gamma_{fs}
\end{align*}
$$


## Vertical boundaries - damping zone
In the model, the waves propagate from left to right. When an incident wave interacts with the floating structure, it will partially reflect at the structure boundary and be partially transmitted. It is assumed that these reflected, and transmitted waves will not further interfere with the structure. They must therefore fully propagate away from the system and it is not allowed for these waves to reflect at the vertical boundaries of the fluid domain. 

Hence, the boundary condition at the inlet of the domain is determined by the predefined expression for the incident wave, $\phi_{\text{inc}}$, where the BC at the outlet is set to zero:

$$
\begin{align*}
    \vec{n}\cdot\nabla\phi &= 0  &&\text{on} \quad \Gamma_{R}\\
    \vec{n}\cdot\nabla\phi &=\vec{n}\cdot\nabla\phi_{\text{inc}} &&\text{on} \quad \Gamma_{L}
\end{align*}
$$

To assure energy dissipation such that these conditions at the vertical boundaries are satisfied, two damping zones (at both the inlet and the outlet) are constructed in the model. To this end, a method by [Kim Woo Min](http://dx.doi.org/10.1016/j.oceaneng.2013.10.012) is used, who introduced two additional terms for the kinematic boundary condition which dissipate the wave energy. This results in the following alteration for the kinematic condition at the free surface:

$$
\begin{equation*}
     -i\omega\eta - \vec{n}\cdot\nabla\phi + μ_1(η - \eta_{\text{inc}} ) + \frac{μ_2}{g}(\phi- \phi_{\text{inc}}) = 0 \qquad \text{on} \quad \Gamma_{fs}
\end{equation*}
$$

Where $\eta_{\text{inc}}$ and $\phi_{\text{inc}}$ (defined by the incident wave expression) are the reference values when the computational domain is not disturbed by any structure during the wave propagation. $μ_1$ and $\mu_2$ are the damping coefficients, which are dependent on each other to ensure that no dispersion occurs:
$$
\begin{equation*}
    \mu_1(x)=\mu_{0}\left[1-\sin \left\{\frac{\pi}{2}\left(\frac{x}{L_{d,in}}\right)\right\}\right] + \mu_{0}\left[1-\cos \left\{\frac{\pi}{2}\left(\frac{x-x_{d}}{L_{d,out}}\right)\right\}\right]
\end{equation*}
$$
$$
\begin{equation*}
    \mu_2(x) = - \frac{\mu_1(x)^2}{4}
\end{equation*}
$$

Where $L_{d,in}$ and $L_{d,out}$ are the total lengths of the respective damping zone. $x_{d}$ is starting point of the damping zone at the oulet. The input value, $\mu_0$, is selected on a trial and error basis, depending on the wave characteristics.


## Fluid-structure interface
Also at the fluid-structure interface both a kinematic and a dynamic condition is defined. For the kinematic condition, the elevation of the water surface is equal to the displacement of VLFS. For the dynamic boundary condition, the water pressure at the surface acts as a distributed load at the bottom of the structure. Hence, the expression for the water pressure is substituted in the EoM of the VLFS. 

Since the equation of motion for the structure is substituted as a complex boundary condition of the fluid domain, it is more convenience to rewrite the vertical displacement of the structure - previously indicated with $w$ - in terms of the surface elevation $\eta$. Therefore, the boundary conditions at the fluid-structure interface yield: 

$$
\begin{align*}
    \text{kinematic:}& \qquad -i\omega\eta  = \vec{n}\cdot\nabla\phi  &&\text{on} \quad \Gamma_{str} \\
    \text{dynamic:}&\qquad -\omega^2\alpha_{b1} \eta + \alpha_{b2}\frac{\partial^4 \eta}{\partial x^4} - i\omega\phi + gη = 0 &&\text{on} \quad \Gamma_{str}
\end{align*}
$$

Where $\alpha_{b1} = \frac{\rho_b h_b}{\rho_w}$ , and $\alpha_{b2} = \frac{EI}{\rho_w}$.


## Structural domain

The moments and shear forces at the free ends of the VLFS should be equal to zero. Consequently, for the first and the last beam element, the following four dynamic BCs must be satisfied:
$$
\begin{align*}
    \frac{\partial^2 \eta}{\partial x^2} = \frac{\partial^3 \eta}{\partial x^3} = 0 \qquad \text{at} \quad x\,=\,0  &&\text{on} \quad \Gamma_{str}\\
    \frac{\partial^2 \eta}{\partial x^2} = \frac{\partial^3\eta}{\partial x^3} = 0 \qquad \text{at} \quad  x\,=\,L &&\text{on} \quad \Gamma_{str}
\end{align*}
$$

Finally, some conditions at the connection points of the VLFS are determined. At each connection, the kinematic conditions (i.e. continuity in displacement and rotation), and the dynamic boundary conditions (i.e. force equilibrium) are satisfied. At each connection location $x_j$ the following 4 ICs should hold: 

$$
\begin{align*}
    \eta \big|_{x_j^-} \, &= \, η \big|_{x_j^+}\\
    \frac{\partial \eta}{\partial x} \big|_{x_j^-} &= \frac{∂η}{\partial x} \big|_{x_j^+}\\
    EI\frac{\partial^2 η}{\partial x^2} \big|_{x_j^-} = EI\frac{\partial^2 η}{\partial x^2}\big|_{x_j^+} &= k_{r}\Bigg( \frac{\partial η}{\partial x}\big|_{x_j^-} - \frac{\partial η}{\partial x} \big|_{x_j^+} \Bigg)\\
    EI\frac{\partial^3 \eta}{\partial x^3} \big|_{x_j^-} &= EI \frac{\partial^3 η}{\partial x^3} \big|_{x_j^+}
\end{align*}
$$


To conclude, the motion of the fluid is described by the Laplace equation within the fluid domain, where the equation of motion for the structure is substituted as a complex boundary condition for the fluid-structure interface. Together with the other boundary conditions as defined above, this results in a well-posed problem, which can now be solved.


## Expressions of the incident wave
According to the linear wave theory, the potential of the incident wave $\phi_{\text{inc}}$ for perpendicular waves and considering water of finite depth, is written in the following form:
$$
\begin{equation*}
    \phi_{\text{inc}} = \frac{gA_w}{i\omega}\,\frac{\cosh{k(z+d)}}{\cosh{kd}}\,\text{e}^{ikx}
\end{equation*}
$$
Where $A_w$ is the wave amplitude of the incident wave. The incident wave number $k$, the wave frequency $\omega$, and the incident wavelength $\lambda$, are related through the dispersion relation:
$$
\begin{equation*}
    gk\tanh(kd) \, = \, \omega^2
\end{equation*}
$$
Where, $k=2\pi/\lambda$ is defined as the wave number


Consequently, the expressions for the surface elevation of the incident wave $\eta_{\text{inc}}$ and the velocity of the fluid in x-direction over the boundary $\Gamma_L$ are formulated as follows:

$$
\begin{align*}
    \frac{\partial\phi_{\text{inc}}}{\partial z} = -i\omega\eta_{\text{inc}} \qquad \rightarrow \qquad  \quad \eta_{\text{inc}} &= \frac{gA_w k}{\omega^2}\,\frac{\sinh{k(z+d)}}{\cosh{kd}}\,\text{e}^{ikx}\\
     \vec{n}\cdot\nabla\phi_{\text{inc}} \big|_{\Gamma_{L}} \qquad \rightarrow \qquad  \frac{\partial\phi_{\text{inc}}}{\partial x} &= \frac{gA_w k}{\omega}\,\frac{\cosh{k(z+d)}}{\cosh{kd}}\,\text{e}^{ikx}
\end{align*}
$$


# Numerical model

Now that the mathematical formulation behind the model has been described, the system of equations can be rewritten into the weak formulation and inserted into the numerical model. However, in order to do so, the input parameters regarding the incident wave and material properties are defined first, afterwhich the numerical model is set up.


## Input parameters

### Material Properties

For this tutorial it is assumed that the VLFS has a solid rectangular cross-section with a height of 2m. The type of connection can be determined by varying the rotational stiffness parameter, where $\xi=0$ corresponds to a hinged connection, and $ξ=650$ to a rigid connection. In this tutorial a value of $ξ=0$ is applied. 

Finally, the length of the separated modules can be altered by changing the locations of the connections with the connection location parameter $β$. Here, it is assumed that the four modules all have the same length, i.e. $β=0.25$.  
 
```julia 
ρ_b = 250               # mass density VLFS [kg m^-3]
h_b = 2                 # height VLFS [m] 
L = 1000                # total length VLFS [m]
I = 1/12*h_b^3          # second moment of inertia per unit meter width (I/b) 
E = 12e9                # Youngs modulus [N/m^2]
EI_b = E*I              # bending stiffness VLFS per unit meter width [Nm/m]  
ξ = 0                   # rotational stiffness parameter
k_r = ξ*EI_b/L          # rotational stiffness connection [Nm]
β = 0.25                # connection location parameter 

α1_b = ρ_b*h_b/ρ_w        
α2_b = EI_b/ρ_w
```

### Fluid domain and incident wave

The fluid characteristics, the dimensions of the fluid domain, and the conditions for the incident wave are described here below; where the input parameters correspond to mild wave conditions.

```julia
g = 9.81                    # gravitational constant
ρ_w = 1025                  # water density [kg m^-3]
d = 30                      # water depth [m]
L_fd = 2*L                  # total length fluid domain [m]

λ = 140                     # incident wave length [m]
k = 2*π / λ                 # wave number 
ω = sqrt(g*k*tanh(k*d))     # dispersion relation                         
A_w = 0.75                  # amplitude of incident wave


ϕ_in(x) = (g*A_w/(1im*ω))*(cosh(k*x[2]) / cosh(k*d))*exp(1im*k*x[1])            # expression flow potential of incident wave
u_in(x) = (g*A_w*k/ω)*(cosh(k*x[2]) / cosh(k*d))*exp(1im*k*x[1])                # expression flow velocity in x-direction over the boundary $\Gamma_L$
η_in(x) = (g*A_w*k/ω^2)*(sinh(k*x[2]) / cosh(k*d))*exp(1im*k*x[1])              # expression of the surface elevation of the incident wave
```

## Domain

Next, the numerical domain is defined using the `Gridap` library, for which the different zones are specified. This includes the inlet damping zone, the free surface zone in front of the VLFS, the VLFS, the free surface zone behind the structure, and the outlet damping zone. For this tutorial, both damping zones, and the VLFS are set equal to L, while the free surface zones are equal to 1/2 L.   

```julia
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces


Ld_in  = 1L                                                 # Inlet damping zone length
Lb     = L                                                  # VLFS length
Ld_out = 1L                                                 # Outlet damping zone length
LΓ = Ld_in+L_fd+Ld_out                                      # length fluid domain inc. damping zones
domain = (0,LΓ,0,d)                                         # fluid domain  

xᵢₙ = Ld_in                                                  # x-coordinate end frontal damping zone
xb₀ = Ld_in + 0.5(L_fd-L)                                    # x-coordinate start of VLFS 
xb₁ = xb₀ + Lb                                               # x-coordinate end of VLFS
xd = xb₁ + 0.5(L_fd-L)                                       # x-coordinate initial outlet damping zone point
xj = [(xb₀ + β*Lb), (xb₀ + 2*β*Lb), (xb₀ + 3*β*Lb)]          # x-coordinates that correspond with locations of the connections
```


### Damping zone 
The functions regarding the damping coefficients are defined with the expressions below, where the initial damping is set to $μ₀ = 10$. 

```julia
μ₀ = 10

μ₁(x::VectorValue) = μ₀* (1.0 - cos(π/2*(x[1]-xd)/Ld_out)) * (x[1] > xd) + μ₀*(1.0 - sin(π/2*x[1]/Ld_in)) * (x[1] < xᵢₙ) 
μ₂(x::VectorValue) = -(μ₁(x)^2)/4

ηd = x -> μ₁(x)*η_in(x)*(x[1]<xᵢₙ)
ϕd = x -> μ₂(x)*ϕ_in(x)*(x[1]<xᵢₙ)                                     
```


### Construct the discrete model of the full domain
To construct the discrete model of the full domain, a partition is defined. The resolution of the horizontal axis should be high enough to obtain an accurate approximation of the propagating wave in x-direction. With respect to the vertical axis, the main interest is focused on the top layer of the fluid domain, near the surface; as the vertical velocity profile converts to zero at the bottom boundary. Therefore, a function is written to force an uneven distribution of the cell spaces in z-direction, such that the resolution is fine at the free surface and becomes coarser towards the bottom of the domain. For this tutorial, $nx=80$ and $nz=20$ are used regarding the partition in x-direction and z-direction, respectively.

Furthermore, to impose the conditions defined at the connections, it is important that the locations of the respective connections, $xj$, coincides with one the nodes of the finite elements. Hence, this is obtained by defining the distribution of the cells as a multiplication of connection location parameter $\beta$.        

Accordingly, the discrete model of the fluid domain $Ω$ is constructed using the function `CartesianDiscreteModel`. The function `simplexify` is used to change the mesh to an affine reference map, which is necessary to have the mapping work.


```julia
nx = 80                               
nz = 20                               
partition = (Int(4/β)*nx,nz)

function f_z(x)                 # function to redistribute refinement of elements in z-direction
  if x == d
      return d
  end
  i = x / (d/nz)
  return d-d/(2^i)
end

map(x) = VectorValue(x[1], f_z(x[2]))
model_Ω=simplexify(CartesianDiscreteModel(domain,partition,map=map))
```

## Boundaries
Now that the discrete fluid domain `model_Ω` has been constructed, the corresponding boundaries of the discrete domain are obtained. First, the initial boundaries of the fluid domain are mapped to their corresponding entity of the fluid domain and labeled using the function `add_tag_from_tags!`.  

```julia
# Define surfaces that can be created without masks (only valid when using rectangular shape and a CartesianDiscreteModel)
labels_Ω = get_face_labeling(model_Ω) 
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)
```

The different zones of the surface boundary are then defined using a `mask`. First, a mask is created for the complete surface boundary of the fluid domain, afterwhich a separated discrete model is created, i.e. `model_Γ`.

```julia
# Masks in Ω
# ==========
# From all the entities of dimension 1 in Ω (edges in the full domain), 
# set to True all the edges that are in Γ and False otherwise.
Γ_mask_in_Ω = get_face_mask(labels_Ω,"surface",1)
# Create a list of indices that have value=True in the mask of active edges in Omega.
Γ_to_Ω_dim1 = findall(Γ_mask_in_Ω)

# Discrete model on the boundary (to construct a FE space only on the boundary Γ)
model_Γ = BoundaryDiscreteModel(Polytope{1},model_Ω,Γ_to_Ω_dim1)
```

Next, auxiliar functions are defined; to check if the elements are either located in one of the damping zones, at the free surface, within the structure, at the boundary of the structure, or at a connection.

```julia 
# Auxiliar functions
# ==================
# Check if an element is inside the beam
function is_beam(coords)
    n = length(coords)
    x = (1/n)*sum(coords)
    (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ d )
end

# Check if an element is inside the inlet damping zone
function is_inlet_damping(xs)
  # xs is a vector of points with the coordinates of a given element of dimension "d"
  n = length(xs) # number of points in each element (in an edge there will be two points)
  x = (1/n)*sum(xs) # take the average (centre of the element)
  (x[1] <= Ld_in ) * ( x[2] ≈ d ) # check if the centre is in the inlet damping zone and on the surface
end

# Check if an element is inside the outlet damping zone
function is_outlet_damping(xs)
  # xs is a vector of points with the coordinates of a given element of dimension "d"
  n = length(xs) # number of points in each element (in an edge there will be two points)
  x = (1/n)*sum(xs) # take the average (centre of the element)
  ( (LΓ - Ld_out) <= x[1] ) * ( x[2] ≈ d ) # check if the centre is in the inlet damping zone and on the surface
end

# Check if an element is on the beam boundary
function is_beam_boundary(xs)
  is_on_xb₀ = [x[1]≈xb₀ for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
  is_on_xb₁ = [x[1]≈xb₁ for x in xs]
  element_on_xb₀ = minimum(is_on_xb₀) # Boolean with "true" if at least one entry is true, "false" otherwise.
  element_on_xb₁ = minimum(is_on_xb₁)
  element_on_xb₀ | element_on_xb₁ # Return "true" if any of the two cases is true
end


# Check if an element is a joint
function is_a_joint(xs)
  is_on_xj = Array{Bool,1}(UndefInitializer(),length(xs))
  is_on_xj .= false
  for xi in xj
    is_on_xj = is_on_xj .| [x[1]≈xi for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
  end
  element_on_xj = minimum(is_on_xj) # Boolean with "true" if at least one entry is true, "false" otherwise.
  element_on_xj 
end
```

The separated zones of the surface boundary are constructed, for which different masks are created with regard to the structure, the connection points, and the free surface. Accordingly, the entities that belong to the respective part of the surface boundary are identified and assigned to a new label.  

```julia 
# Masks in Γ
# ==========
labels_Γ = get_face_labeling(model_Γ) # get the face labeling of model_Γ
topo = get_grid_topology(model_Γ) # get the topology of model_Γ (the information about how the model is defined)
D = num_cell_dims(model_Γ) # spatial dimension of the model_Γ

# Construct the mask for the beam
entity_tag_beam = num_entities(labels_Γ) + 1 # create a new tag for the beam
for d in 0:D # loop over dimensions
  grid_Γ_dim_d = Grid(ReferenceFE{d},model_Γ) # construct a grid from the entities of dimension "d" in Γ
  coords_grid_Γ_dim_d = get_cell_coordinates(grid_Γ_dim_d) # get the coordinates of entities of dimension "d" of the grid (array of vectors of points). If model_Γ has 10 edges, this will be an array of 10 entries where each entry will contain a vector of 2 points.
  beam_mask_in_grid_Γ_dim_d = lazy_map(is_beam,coords_grid_Γ_dim_d) # beam mask with the entities of dimension "d" 
  beam_boundary_mask_in_grid_Γ_dim_d = lazy_map(is_beam_boundary,coords_grid_Γ_dim_d) # beam boundary mask with the entities of dimension "d"
  joint_mask_in_grid_Γ_dim_d = lazy_map(is_a_joint,coords_grid_Γ_dim_d) # joint mask with the entities of dimension "d"
  # Create a list of indices that have value=True in the mask of active entities of dimension "d" in Gamma. Excluding the entities on the beam boundary and on the joint.
  beam_to_Γ_dim_d = findall( beam_mask_in_grid_Γ_dim_d .&
                             .!beam_boundary_mask_in_grid_Γ_dim_d .&
                             .!joint_mask_in_grid_Γ_dim_d)
  # Tag the faces of dimension "d" (all the faces indexed in beam_to_Γ_dim_d)
  for face in beam_to_Γ_dim_d
    labels_Γ.d_to_dface_to_entity[d+1][face] = entity_tag_beam
  end
end

# Add a name to the tag in labels_Γ
add_tag!(labels_Γ,"beam",[entity_tag_beam])


# Construct the mask for the joint
entity_tag_joint = entity_tag_beam + 1
# Here we don't loop because we'll only have entities of dimension 0 (points)
grid_Γ_dim_0 = Grid(ReferenceFE{0},model_Γ)
coords_grid_Γ_dim_0 = get_cell_coordinates(grid_Γ_dim_0)
joint_mask_in_grid_Γ_dim_0 = lazy_map(is_a_joint,coords_grid_Γ_dim_0)
# Create a list of indices of points in Γ that are joint
joint_to_Γ_dim_0 = findall( joint_mask_in_grid_Γ_dim_0)
# Tag the points (all the points indexed in joint_to_Γ_dim_0)
for point in joint_to_Γ_dim_0
  labels_Γ.d_to_dface_to_entity[1][point] = entity_tag_joint
end
# Add a name to the tag in labels_Γ
add_tag!(labels_Γ,"joint",[entity_tag_joint])
```

Next, the masks of edges in `Ω`, and the masks of points in `Γ` are obtained, which will be used to construct the triangulations for the floating structure, the free surface, and the set of internal points of the structure (skeleton).

```julia
# Mask of edges only in Ω
Γ_mask_in_Ω_dim_1 = get_face_mask(labels_Ω,"surface",1)         # get the mask of edges in Ω that are on Γ
Γ_to_Ω_dim_1 = findall(Γ_mask_in_Ω_dim_1)                       # Get indices of edges of Ω that are on Γ
Γstr_mask_in_Γ_dim_1 = get_face_mask(labels_Γ,"beam",1)         # get the mask of edges in Γ that are on the beam
Γstr_to_Γ_dim_1 = findall(Γstr_mask_in_Γ_dim_1)                 # get indices of edges of Γ that are in the beam
Γfs_to_Γ_dim_1 = findall(!,Γstr_mask_in_Γ_dim_1)                # get indices of edges of Γ that are in the free surface
Γstr_to_Ω_dim_1 = view(Γ_to_Ω_dim_1,Γstr_to_Γ_dim_1)            # get indices of edges of Ω that are in the beam. To create that we "concatenate" two sets (from Γstr to Γ and from Γ to Ω)
Γfs_to_Ω_dim_1 = view(Γ_to_Ω_dim_1,Γfs_to_Γ_dim_1)              # Idem for the free surface

# Mask of points only in Γ
Λb_mask_in_Γ_dim_0 = get_face_mask(labels_Γ,"beam",0)
Λj_mask_in_Γ_dim_0 = get_face_mask(labels_Γ,"joint",0)
```


## Triangulations
Accordingly, the triangulations for the domain, the boundaries, and the interior points of the structure can be easily obtained

```julia
Ω = Triangulation(model_Ω) # triangulation of the full domain
Γ = Triangulation(model_Γ) # triangulation of the boundary (free_surface + beam)

# Now we can construct sub-triangulations
Γstr = BoundaryTriangulation(model_Ω,Γstr_to_Ω_dim_1)
Γfs = BoundaryTriangulation(model_Ω,Γfs_to_Ω_dim_1)
Γin = BoundaryTriangulation(model_Ω, tags = ["inlet"])

# Now construct a skeleton triangulation for the beam and joint
Λb = SkeletonTriangulation(model_Γ,Λb_mask_in_Γ_dim_0)
Λj = SkeletonTriangulation(model_Γ,Λj_mask_in_Γ_dim_0)
```



## Quadratures 
Finally, the quadratures are specified for the domain, the boundaries, and interior points of the structure

```julia
order = 2
dΩ = Measure(Ω,2*order)
dΓ_str = Measure(Γstr,2*order)
dΓ_fs = Measure(Γfs,2*order)
dΓ_in = Measure(Γin,2*order)
nΛb = get_normal_vector(Λb)
dΛb = Measure(Λb,2*order)
nΛj = get_normal_vector(Λj)
dΛj = Measure(Λj,2*order)
```

## FE spaces
As the numerical domain, and the specific boundaries have been defined, the test spaces can be
constructed. To this end, two spaces are built with regard to the internal domain `model_Ω` and the surface `model_Γ`, for which linear lagrangian shape functions are used as reference element. Subsequently, the trial spaces are obtained from the test spaces. The separated spaces are then combined for the full numerical domain using the function `MultiFieldFESpace`.

```julia
reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(model_Ω,reffe, vector_type=Vector{ComplexF64}  )
  V_Γ = TestFESpace(model_Γ,reffe, vector_type=Vector{ComplexF64}  )
  U_Ω = TrialFESpace(V_Ω)
  U_Γ = TrialFESpace(V_Γ)
  Y = MultiFieldFESpace([V_Ω,V_Γ])
  X = MultiFieldFESpace([U_Ω,U_Γ])
```



# Weak form & CG/DG approach

The final step to solve the FSI problem with Gridap is to obtain the bilinear form of the mathematical formulation, $a((ϕ,η),(w,v))$ (see below). The bilinear form contains both first, and second order derivatives. The finite element spaces therefore require continuous gradients between elements, i.e. C1 continuity across elements. As proposed by Colomés et al., this can be achieved by applying the Continuous Galerkin / Discontinuous Galerkin (CG/DG) approach for fourth order operators, where the discrete functions are continuous at the element nodes, but the gradient is discontinuous. Therefore, linear Lagrangian elements can be applied, while continuity of the gradient over adjacent elements is weakly enforced by means of an interior penalty approach. 

To obtain the bilinear form, the following steps are applied:

- The first row is obtained by multiplying the Laplace equation by the test function, $w$, after which it is integrated over the domain $Ω$ and then integrated by parts. 
- The normal velocity at the boundaries is replaced by the respective kinematic boundary condition of the free surface $\Gamma_{fs}$, see second row, where the conditions at the seabed $\Gamma_b$, and the vertical boundary $\Gamma_R$ are canceled out as they go to zero.
- The third row includes the dynamic boundary condition on the free surface. According to the monolithic approach described by [Akkerman et al.](https://doi.org/10.1016/j.oceaneng.2020.107114), here the dynamic condition is multiplied by a modified test function, $v+\alpha_{fs}w$, and integrating over the free surface boundary $\Gamma_{fs}$, where the term $\alpha_{fs}w$ is added to the test function $v$ to guarantee coercivity of the system. 
- The fourth row includes the dynamic boundary condition on the beam surface, multiplied by test function $v$, integrated over the fluid-structure interface $\Gamma_{str}$. In this case, however, the fourth order term is integrated by parts twice, where the resulting integrals over the structure boundaries are canceled out again, as it was assumed that the conditions at the free ends of the VLFS are equal to zero.
- Implementing the interior penalty approach in a similar way as discussed in the Gridap tutorial [Poisson equation (with DG)](https://gridap.github.io/Tutorials/stable/pages/t006_dg_discretization/), the third, and second last row weakly forces continuity of the surface elevation gradients between structural elements. 
- Ultimately, the approach to weakly force conditions between adjacent elements is also used to impose the interface condition regarding the bending moment at the connections, as can be seen in the final row. 



## weak form 
Hence, the complete bilinear form is described below, where the linear form involves the boundary conditions which are defined by the known expressions for $\phi_{\text{inc}}$, and $\eta_{\text{inc}}$.

```julia
const h = β*Lb/nx
const γ_m = 1.0e2*order*(order+1)
const αh = -1im
const βh_fs = 0.5
αh_fs = αh*ω/g*(1-βh_fs)/βh_fs



a((ϕ,η),(w,v)) =      ∫(  ∇(w)⋅∇(ϕ) )dΩ   +   
                      ∫(  (1im*ω*w*η)  - μ₁*η*w - μ₂*ϕ*w/g )dΓ_fs   +
                      ∫(  βh_fs*(v + αh_fs*w)*g*η  +   βh_fs*(-1im*ω)*(v + αh_fs*w)*ϕ )dΓ_fs   +
                      ∫(  (-ω^2*α1_b + g)*v*η +  Δ(v)*(α2_b*Δ(η)) +  (-1im*ω*v*ϕ)   +   (1im*ω*w*η)  )dΓ_str   +
                      ∫(   - (jump(∇(v)⋅nΛb) * mean(α2_b*Δ(η))) - (mean(Δ(v)) * jump(α2_b*∇(η)⋅nΛb))  + 
                          γ_m/h*( jump(∇(v)⋅nΛb) * jump(α2_b*∇(η)⋅nΛb))  )dΛb +
                      ∫(  (1/ρ_w)*(jump(∇(v)⋅nΛj) * k_r * jump(∇(η)⋅nΛj)) )dΛj
                      

l((w,v)) =            ∫( w*u_in )dΓ_in - ∫( w*ηd + w*ϕd/g )dΓ_fs
```


# Solver
Ultimately, with the function `AffineFEOperator`, the equations are assembled into a matrix, where the numerical solver `Gridap.solve` is used to find the approximated solution for `ϕh`, and `ηh`. 

The bending moment of the VLFS is defined by the second derivative of the displacement with respect to x, multiplied by the bending stiffness EI. Accordingly, the moment distribution along the structure can be obtained from the solution for `ηh`; using  `interpolate_everywhere`  over the surface boundary, `V_Γ`.  


```julia
op = AffineFEOperator(a, l, X, Y)

ϕh , ηh = Gridap.solve(op)
Moment = interpolate_everywhere(EI_b*Δ(ηh),V_Γ)
```

# Visualisation
The results can be inspected by writing it into a `vtk file`, which will generate a file that contains the real and imaginay part of the problem solutions.

```julia
writevtk(Ω, "FSI_multiVLFS_phi", cellfields=["phi_re"=>real(ϕh), "phi_im"=>imag(ϕh)])
writevtk(Γ, "FSI_multiVLFS_eta", cellfields=["eta_re"=>real(ηh), "eta_im"=>imag(ηh)])
writevtk(Γ, "FSI_multiVLFS_moment_distribution", cellfields=["M_re"=>real(Moment), "M_im"=>imag(Moment)])
```





