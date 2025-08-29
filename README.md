# `wolfd2`: 2D Incompressible Flow

`wolfd2` simulates 2D incompressible flow by numerically solving the
Navier-Stokes equations over generalized coordinates. Pressure-velocity coupling
is enforced using Gresho's 'Projection 1' method. Spatial derivatives are
discretized over finite volumes using cell-based indexing.  Time derivatives 
are integrated using Douglas-Gunn ADI time-splitting.  

The program also solves a thermal energy equation.  Buoyancy source terms are
included in the momentum equations. The modified Shuman filter is used to
control aliasing and cell-Reynolds (or cell-Peclet) number oscillations.

Turbulet flow is simulated using the Additive Turbulent Decomposition (ATD)
method of Hylin and McDonough.  This implementation of ATD uses linear
combinations of (appropriately-scaled) logistic maps to simulate the small-scale
velocity fluctuations.

The program simulates porous media flow by solving the momentum and thermal
energy equations with Darcy and Forschheimer terms.

Particle trajectories are calculated using a Lagrangiatn approach.

This program was written as part of masters thesis research under the
supervision of Prof. James M. McDonough, Computational Fluid Dynamics
Laboratory, Dept. of Mechanical Engineering, University of Kentucky.

Version 0.3

June 13, 1999
