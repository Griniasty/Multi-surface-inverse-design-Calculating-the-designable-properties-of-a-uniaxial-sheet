# Multi-valued-inverse-design-multiple-3D-surfaces-from-one-flat-sheet
A single uniaxial sheet, with N designbale features, can be inverse designed to morph into N+1 target surface geometries.
The system of partial differential equations for solving this inverse problem is soon to be published [1].
The Mathematica package, inver and two notebooks in this repository implement an Euler type numerical solution of this inverse problem.
The inverse design 
The first notebook, "Fig1.nb", reproduces Figure 1. of [1], where a single sheet morphs into two surface geometries.
In it a warp and weft integration scheme is implemented, allowing only for analytic specification of the surface geometries.
The second notebook, "Fig3.nb", reproduces Figure 3. of [1], where a single sheet morphs into three surface geometries, in response to two independent stimuli.
In it a diagonal integration scheme is implemeted, allowing for both analytic and numerical specification of the surface geometries.

[1] Itay Griniasty, Cyrus Mostajeran and Itai Cohen, "Multi-valued inverse design: multiple 3D surfaces from one flat sheet", to be published.
