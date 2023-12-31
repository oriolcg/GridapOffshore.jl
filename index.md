@def title = "**Gridap Offshore**: A library of Finite Element tutorials for/by Offshore and Hydraulic Engineers"
@def tags = ["syntax", "code"]
@def mintoclevel=1
@def maxtoclevel=1

This is a library of tutorials made for and by Hydraulic and Offshore Engineering students at [**Delft University of Technology**](www.tudelft.nl). The goal of this library is two fold: 

1. Guide students on the use of [*Gridap.jl*](https://github.com/gridap/Gridap.jl) to solve PDEs with the Finite Element method, providing tutorials for a wide variety of background levels.
2. Give visibility to the work of MSc students that used Gridap.jl in their MSc thesis/project.

The collection of tutorials posted in this site are based on [*Gridap.jl*](https://gridap.github.io/Gridap.jl/stable/) a pure *Julia* Finite Element library. Here you will find tutorials covering contents in computational modeling taught in the [Civil Engineering MSc at TU Delft](https://www.tudelft.nl/en/education/programmes/masters/civil-engineering/msc-civil-engineering), and tutorials related to MSc thesis from different masters, including:
- [Master of Offshore and Dredging Engineering](https://www.tudelft.nl/onderwijs/opleidingen/masters/offshore-dredging-engineering/msc-offshore-dredging-engineering)
- [Master of Civil Engineering: Hydraulic Engineering](https://www.tudelft.nl/onderwijs/opleidingen/masters/ce/msc-civil-engineering-test/old/old/old/oud/old/old/old/msc-programme/track-hydraulic-engineering)

List of tutorials:

\toc

# Theory tutorials

# MSc Thesis tutorials

- January 2022: [~~~<b>Waves through porous medium</b>~~~](/Fluids/Porous/porous) by Joël Ruesen

This tutorial shows how wave-progression through a porous medium is modelled. The model uses viscous incompressible Navier Stokes in combination with Darcy-Forchheimer resistance terms in the momentum balance, implemented using the Gridap library.

~~~<u>Reference</u>~~~: Ruesen, Joël. *Wave damping by large-scale offshore kelp farms - a numerical modelling framework using a porous medium approach.* (2022). [MSc thesis](http://resolver.tudelft.nl/)

***
    
- October 2021: [~~~<b>Very Large Floating Structures (VLFS)</b>~~~](/FSI/VLFS/VLFS) by Dorette Regout

This tutorial shows how to solve a Fluid Structure Interaction (FSI) problem using Gridap and provides the instructions to build a 2D model considering a multi-module VLFS, solved in the frequency domain.
@@im
![](/FSI/VLFS/img/numerical_model.png) 
@@
~~~<u>Reference</u>~~~: Regout, Dorette. *Hydroelastic Analysis of a Multi-Module Very Large Floating Structure.* (2021). [MSc thesis](http://resolver.tudelft.nl/uuid:838e2e14-d0e8-49dc-bad4-e7e132b248bc)

***

- July 2021: [~~~<b>Very Flexible Floating Structures (VFFS)</b>~~~](/FSI/VFFS/VFFS) by Sjoerd van Hoof

This tutorial shows how a Fluid Structure Interaction (FSI) in a 2D domain is modelled. Potential flow is used to model the fluid and on top a Dynamic Euler-Bernoulli beam is located that serves as the floating structure.
@@im-50
![](/FSI/VFFS/img/viridis_3D.png) 
@@    
~~~<u>Reference</u>~~~: van Hoof, Sjoerd. *Hydroelastic wave deformation of Very Flexible Floating Structures: A performance study of a monolithic finite element model.* (2021). [MSc thesis](http://resolver.tudelft.nl/uuid:6652a9ee-61a6-4c4a-9b28-ef4fda1010f9)