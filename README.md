# unsat_flow

**UG4-App** implementing a density driven flow equation system

## Documentation
The equations used are: </br>
<img src="https://render.githubusercontent.com/render/math?math=\partial_t (\Phi \rho_w S_w) + \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w"> </br>
<sub>($\partial_t (\Phi \rho_w S_w) + \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$)</sub> </br>
and </br>
<img src="https://render.githubusercontent.com/render/math?math=\partial_t (\Phi \rho_w S_w \omega) + \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] = \rho_w \omega \Gamma_w"> </br>
<sub>($\partial_t (\Phi \rho_w S_w \omega) + \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] = \rho_w \omega \Gamma_w$)</sub>

The system is implemented using UG4's luashell, the Richards and ConvectionDiffusion plugins.

Usage:


## Examples

## Dependencies
This app depends on the Richards Plugin by Dr. Arne Nägel and UG4's ConvectionDiffusion Plugin.

## References
[1] Eckhard Schneid: Hybrid-Gemischte Finite-Elemente-Diskretisierung der Richards-Gleichung, Dissertation, FAU Erlangen, 2000

[2] K. Johannsen: Numerische Aspekte dichtegetriebener Strömung in porösen Medien, Habilitationsschrift, Universität Heidelberg, 2004

[3] P. Frolkovic: Application of level set method for groundwater flow with moving boundary. Advances in Water Resources 47 (2012) 56–66
