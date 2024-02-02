@def title = "<b>Gridap Offshore</b><br>A library of Finite Element tutorials for/by Offshore and Hydraulic Engineers"
@def tags = ["syntax", "code"]
@def mintoclevel=1
@def maxtoclevel=1

This is a library of tutorials made for and by Hydraulic and Offshore Engineering students at [**Delft University of Technology**](www.tudelft.nl). The goal of this library is two fold: 

1. Guide students on the use of [*Gridap.jl*](https://github.com/gridap/Gridap.jl) to solve PDEs, providing tutorials for a wide variety of background levels.
2. Give visibility to the work of MSc students that used Gridap.jl in their MSc thesis/project.

The collection of tutorials posted in this site are based on [*Gridap.jl*](https://gridap.github.io/Gridap.jl/stable/) a pure *Julia* Finite Element library. Here you will find tutorials covering contents in computational modeling taught in the [Civil Engineering MSc at TU Delft](https://www.tudelft.nl/en/education/programmes/masters/civil-engineering/msc-civil-engineering), and tutorials related to MSc thesis from different masters, including:
- [Master of Offshore and Dredging Engineering](https://www.tudelft.nl/onderwijs/opleidingen/masters/offshore-dredging-engineering/msc-offshore-dredging-engineering)
- [Master of Civil Engineering: Hydraulic Engineering](https://www.tudelft.nl/onderwijs/opleidingen/masters/ce/msc-civil-engineering-test/old/old/old/oud/old/old/old/msc-programme/track-hydraulic-engineering)

## List of tutorials:

\toc

# Theory tutorials

1. [~~~<b>Solving PDEs with Gridap.jl</b>~~~](./Theory/tutorial_1/intro_FE_1D)
1. [~~~<b>Solving the Euler-Bernoulli equation with Continuous/Discontinuos FEs</b>~~~](./Theory/tutorial_EulerBernoulli/euler_bernoulli)
1. [~~~<b>Solving the Timoshenko beam equation: approaches to avoid shear locking</b>~~~](./Theory/tutorial_Timoshenko/Timoshenko)
1. [~~~<b>Strong vs weak Dirichlet boundary conditions</b>~~~](./Theory/tutorial_Poisson_weakBCs/Poisson_weakBC)

# MSc Thesis tutorials

## [~~~<b>Waves through porous medium</b>~~~](./MSc_thesis/Porous/2022_J_Ruesen/porous) by Joël Ruesen, January 2022

This tutorial shows how wave-progression through a porous medium is modelled. The model uses viscous incompressible Navier Stokes in combination with Darcy-Forchheimer resistance terms in the momentum balance, implemented using the Gridap library.

~~~<u>Reference</u>~~~: Ruesen, Joël. *Wave damping by large-scale offshore kelp farms - a numerical modelling framework using a porous medium approach.* (2022). [MSc thesis](http://resolver.tudelft.nl/uuid:fcba6da4-5d83-415d-a5b9-28fc054e7b15)

## [~~~<b>Very Large Floating Structures (VLFS)</b>~~~](./MSc_thesis/Floating/2021_D_Regout/VLFS) by Dorette Regout, October 2021

This tutorial shows how to solve a Fluid Structure Interaction (FSI) problem using Gridap and provides the instructions to build a 2D model considering a multi-module VLFS, solved in the frequency domain.
@@im
![](/MSc_thesis/Floating/2021_D_Regout/img/numerical_model.png) 
@@
~~~<u>Reference</u>~~~: Regout, Dorette. *Hydroelastic Analysis of a Multi-Module Very Large Floating Structure.* (2021). [MSc thesis](http://resolver.tudelft.nl/uuid:838e2e14-d0e8-49dc-bad4-e7e132b248bc)

## [~~~<b>Very Flexible Floating Structures (VFFS)</b>~~~](./MSc_thesis/Floating/2021_S_van_Hoof/VFFS) by Sjoerd van Hoof, July 2021

This tutorial shows how a Fluid Structure Interaction (FSI) in a 2D domain is modelled. Potential flow is used to model the fluid and on top a Dynamic Euler-Bernoulli beam is located that serves as the floating structure.
@@im-50
![](/MSc_thesis/Floating/2021_S_van_Hoof/img/viridis_3D.png) 
@@    
~~~<u>Reference</u>~~~: van Hoof, Sjoerd. *Hydroelastic wave deformation of Very Flexible Floating Structures: A performance study of a monolithic finite element model.* (2021). [MSc thesis](http://resolver.tudelft.nl/uuid:6652a9ee-61a6-4c4a-9b28-ef4fda1010f9)

# Additional material

The tutorials in this library are not covering all the features of Gridap. If you are interested in additional Gridap-related material, please take a look at the following sources:

## Sources to learn Julia: 

- Download the [binary and documentation](https://julialang.org/) 
- [Short free courses](https://juliaacademy.com/courses?preview=logged_out)
- Introduction workshop to julia: [Introduction to Julia | Jose Storopoli | JuliaCon 2022](https://www.youtube.com/watch?v=uiQpwMQZBTA) 

## Sources about Gridap.jl: 

- [Gridap.jl](https://github.com/gridap/Gridap.jl): a FE package to solve PDEs in Julia 
- [Gridap wiki](https://github.com/gridap/Gridap.jl/wiki). In this wiki-page you will find interesting descriptions on: 
  - Getting started with Gridap 
  - Using VS-Code as Julia IDE 
  - REPL workflows 
  - "How-to"  tutorials on creating a julia package, generating documentation, precompiling a package and registering a package in the Julia server. 
- [Gridap tutorials](https://gridap.github.io/Gridap.jl/dev/)

## Presentations at JuliaCon: 

- [Solving partial differential equations in Julia with Gridap.jl | Francesc Verdugo | JuliaCon 2020](https://www.youtube.com/watch?v=txcb3ROQBS4)
- [New tools to solve PDEs in Julia with Gridap.jl | Francesc Verdugo et al | JuliaCon2021](https://www.youtube.com/watch?v=hsQiFP4S5RY)
- [Solving transient PDEs in Julia with Gridap.jl | Oriol Colomes | JuliaCon 2022](https://www.youtube.com/watch?v=heeiSoKnlUk) 
