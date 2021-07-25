# Multi-valued-inverse-design-multiple-3D-surfaces-from-one-flat-sheet
A single uniaxial sheet, with N designable features, can be inverse designed to morph into N+1 target surface geometries.
The system of partial differential equations for solving this inverse problem is published in [1].
The Mathematica package, "uniaxialintegratorV2.wl", contains definitions for functions that numerically integrate these equations.
The code can be generically applied to inverse design any uniaxial sheet. However, it is effective in producing a solution with substantial coverage only when each designable feature, \epsilon, affects only the deformation magnitude along or across the deformation's anisotropic axis, denoted in the paper by the director, n.

Example implementations of the code, where a specific inverse problem is solved, are given in the notebooks "fig1reproduction.nb" and "fig3reproduction.nb"
The first notebook, "fig1reproductionV2.nb", reproduces a semblence of Figure (1) of [1], where a single sheet morphs into two surface geometries.
The data given in this notebook is fully analytic and so it's integration is swift.
The second notebook, "fig3reproductionV2.nb", reproduces Figure (3) of [1], where a single sheet morphs into three surface geometries, in response to two independent stimuli. 

In both notebooks a diagonal integration scheme is implemented, allowing for both analytic and numerical specification of the surface geometries.

[1] Itay Griniasty, Cyrus Mostajeran and Itai Cohen, "Multi-valued inverse design: multiple surface geometries from one flat sheet", arXiv:2102.07840 

A comment:
The current implementation of the digonal integration scheme does not maximize the integrable (u,v) domain; 
The underlying PDEs imply that the boundary is block diagonal in the (u,v) domain. 
The current implementation realizes instead diagonal boundaries.
The solution's coverage can then naturally be improved on by obtaining the theoretically predicted boundary of the domain.
The boundary may be optimized in future versions by introduction of varying step sizes and mesh refinement at singularities.
A perliminerary demonstration is shown in version V1 in the file "fig1reproductionV2.nb", where a warp and weft integration scheme is implemented.
