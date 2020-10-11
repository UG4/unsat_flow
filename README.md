# unsat_flow

**UG4-App** implementing an unsaturated density driven flow equation system

## Documentation
The equations used are: </br>
<img src="https://render.githubusercontent.com/render/math?math=\partial_t (\Phi \rho_w S_w) %2B \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w"> </br>
<sub>($\partial_t (\Phi \rho_w S_w) + \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$)</sub> </br>
and </br>
<img src="https://render.githubusercontent.com/render/math?math=\partial_t (\Phi \rho_w S_w \omega) %2B \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] = \rho_w \omega \Gamma_w"> </br>
<sub>($\partial_t (\Phi \rho_w S_w \omega) + \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] = \rho_w \omega \Gamma_w$)</sub>

Usage:
ugshell -ex unsat_flow_app/unsat_flow_driver.lua --problem-id "example" --check

## Possible arguments
* --problem-id specifies the problems config file, standard: "trench2D"
* --numPreRefs specifies the number of refinements before distribution
* --numRefs specifies the number of refinements after distribution
* --check checks if a problem file has the correct layout

## Examples

## Dependencies
This app depends on UG4's ConvectionDiffusion Plugin.

## References
[1] Eckhard Schneid: Hybrid-Gemischte Finite-Elemente-Diskretisierung der Richards-Gleichung, Dissertation, FAU Erlangen, 2000

[2] K. Johannsen: Numerische Aspekte dichtegetriebener Strömung in porösen Medien, Habilitationsschrift, Universität Heidelberg, 2004

[3] P. Frolkovic: Application of level set method for groundwater flow with moving boundary. Advances in Water Resources 47 (2012) 56–66
