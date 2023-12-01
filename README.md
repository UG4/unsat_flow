# Flow and Transport in unsaturated media

**UG4-App** implementing an unsaturated density driven flow equation system.

Authors: Niklas Conen, Arne Nägel

https://github.com/Nordegraf/unsat_flow


# Model
The Model consists of two coupled PDEs:

$\frac{\partial}{\partial t} (\phi S_w \rho_w) + \nabla \rho_w \vec{q} = \rho_w \Gamma$
$\frac{\partial}{\partial t}(\phi S_w \rho_w \omega) + \nabla  (\vec{q} \rho_w \omega - \phi S_w \rho_w D \nabla \omega) = \phi S_w \rho_w \Gamma $


## Documentation
The equations used are: 

$\partial_t (\Phi \rho_w S_w) + \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$
$\partial_t (\Phi \rho_w S_w \omega) + \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] = \rho_w \omega \Gamma_w$


The first equation models groundwater flow in a porous medium, the second equation models advection and diffusion of a contaminant.

Saturation and hydraulic conductivity are calculated using the van Genuchten Model.

Usage:

ugshell -ex unsat_flow_app/unsat_flow_driver.lua --problem-id "example" --check

## Dependencies
This app depends on the following UG4 plugins: ConvectionDiffusion, LIMEX, Richards.



# Setup
Please refer to the offical documentation of [`ughub`](https://github.com/UG4/ughub) on how to setup and compile UG4 correctly.
This app depends on UG4's ConvectionDiffusion and Richards Plugins. Please compile them with UG4.


Once compiled just place this repository in the apps directory of your UG4 directory. 

# Usage
The app is split into a driver and an util script. Problems are described using configuration scripts. Multiple examples are given in the configs folder. 
To start a simulation the driver needs to be executed using the ugshell:

ugshell -ex unsat_flow_app/unsat_flow_driver.lua --problem-id "example"

## Arguments
Several additional arguments can be set:
* --problem-id: path to the problems config file
* --numPreRefs: number of refinements before distribution
* --numRefs: number of refinements after distribution
* --dt: minimal timestep size
* --newton: solves the problem using a newton solver if set


# References
[1] Bear, J., and Cheng, A. H.-D. Modeling Groundwater Flow and Contaminant Transport. Theory and Applications of Transport in Porous Media. Springer Science+Business Media B.V, 2010.

[2] Van Genuchten, M. A Closed-form Equation for Predicting the Hydraulic Conductivity of Unsaturated Soils. Soil Science Society of America Journal 44 (09 1980).

[3] Richards, L. A. Capillary conduction of liquids through porous mediums. Physics 1, 5 (1931), 318–333.

[4] P. Frolkovic: Application of level set method for groundwater flow with moving boundary. Advances in Water Resources 47 (2012) 56–66

[5] N. Conen: Hydrochemische Modellierung des Stickstoffeintrags durch landwirtschaftliche Nutzflächen in Regionen mit Trinkwasserbrunnen. B.Sc. thesis, Universtität Frankfurt, 2022.

