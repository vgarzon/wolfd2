# `wolfd2`: 2D incompressible flow solver

`wolfd2` solves the incompressible Navier-Stokes equations in generalized
coordinate form. The pressure-velocity coupling is treated via Gresho's
'Projection 1' method. The spatial derivatives are discretized on finite volumes
using a cell-based indexing. The integration of time derivatives is carried out
by applying the Douglas-Gunn ADI time-splitting method.  In addition, the
thermal energy equation can be solved and buoyancy terms can be included as
source terms in the momentum equations. The modified Shuman filter is used to
control aliasing and cell-Reynolds (or cell-Peclet) number oscillations.

Turbulet flow is simulated following the Additive Turbulent Decomposition (ATD)
method of Hylin and McDonough to model turbulent fluctuations. The version of
ATD implemented herein uses linear combinations of (properly-scaled) chaotic
maps to simulate the small-scale part of the solution.

The program simulates porous media flow by solving the momentum and thermal
energy equations with Darcy and Forschheimer terms.

Particle trajectories are calculated using a Lagrangiatn approach.

This program was written as part of masters thesis research under the
supervision of Prof. James M. McDonough, Computational Fluid Dynamics
Laboratory, Dept. of Mechanical Engineering, University of Kentucky

Version 0.3

June 13, 1999
