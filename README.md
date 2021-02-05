# Multi-valued-inverse-design-multiple-3D-surfaces-from-one-flat-sheet
A single uniaxial sheet, with N designable features, can be inverse designed to morph into N+1 target surface geometries.
The system of partial differential equations for solving this inverse problem is soon to be published [1].
The Mathematica package, "uniaxialinversedesignv1.wl", contains definitions for functions that numerically integrate these equations.
The code can be generically applied to inverse design any uniaxial sheet. However, it is effective in producing a solution with substantial coverage only when each designable feature, \epsilon, affects only the deformation magnitude along or across the deformations anisotropic axis, denoted in the paper by the director, n.

Example implementations of the code, where a specific inverse problem is solved, are given in the notebooks "fig1reproduction.nb" and "fig3reproduction.nb"
The first notebook, "fig1reproduction.nb", reproduces Figure 1. of [1], where a single sheet morphs into two surface geometries.
In it a warp and weft integration scheme is implemented, allowing only for analytic specification of the surface geometries.
The second notebook, "fig3reproduction.nb", reproduces Figure 3. of [1], where a single sheet morphs into three surface geometries, in response to two independent stimuli. In it a diagonal integration scheme is implemented, allowing for both analytic and numerical specification of the surface geometries.

[1] Itay Griniasty, Cyrus Mostajeran and Itai Cohen, "Multi-valued inverse design: multiple 3D surfaces from one flat sheet", to be published.
