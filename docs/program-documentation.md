# wolfd2 Program Documentation

## Overview

`wolfd2` is a 2D incompressible flow simulation program that solves the Navier-Stokes equations using advanced numerical methods. The program was developed as part of master's thesis research at the Computational Fluid Dynamics Laboratory, University of Kentucky.

### Key Features

- **Gresho's Projection 1 Method**: Pressure-velocity coupling for incompressible flow
- **Finite Volume Discretization**: Cell-based indexing on generalized coordinates
- **Douglas-Gunn ADI Time-Splitting**: Efficient time integration method
- **Additive Turbulent Decomposition (ATD)**: Turbulence modeling using chaotic maps
- **Lagrangian Particle Tracking**: Particle trajectory calculations
- **Porous Media Support**: Flow through porous media modeling
- **Thermal Energy Equation**: Optional temperature field solution

## Program Architecture

### Main Components

The program consists of the following main modules:

1. **main.f** - Main program entry point and time-stepping loop
2. **parse.f** - Input parameter parsing and file reading
3. **grid.f** - Grid reading and metric computation
4. **bound_cond.f** - Boundary condition setup and application
5. **momentum.f** - Momentum equation solving (auxiliary velocities)
6. **pressure.f** - Pressure Poisson equation solver
7. **thermal.f** - Thermal energy equation solver
8. **small_scale.f** - ATD small-scale turbulence modeling
9. **traject.f** - Lagrangian particle tracking
10. **file_manip.f** - File I/O operations
11. **utility.f** - Utility functions (filtering, averaging, projection)
12. **string.f** - String manipulation functions

## Execution Flow

### Initialization Phase

1. **Parameter Reading** (`ReadParam` in `parse.f`)
   - Reads `input.dat` file
   - Parses simulation parameters (time steps, tolerances, flags)
   - Reads boundary condition file name
   - Sets up output file prefixes

2. **Grid Setup** (`Grid` in `grid.f`)
   - Reads grid file (Plot3D format)
   - Computes metric tensor components
   - Calculates Jacobian determinants
   - Sets up staggered grid coordinates

3. **Boundary Condition Setup** (`SetUpBCs` in `bound_cond.f`)
   - Reads boundary condition file
   - Sets up region types (internal, porous, blockage)
   - Initializes boundary condition values
   - Configures thermal regions

4. **Initial Conditions** (`InitialCond` in `bound_cond.f`)
   - Cold start: Initialize all fields to zero or specified values
   - Restart: Read previous solution from restart file

5. **Initial Mass Conservation** (if cold start)
   - Apply velocity boundary conditions
   - Solve pressure Poisson equation
   - Project velocities to enforce divergence-free condition

6. **Small-Scale Initialization** (if ATD enabled)
   - Initialize small-scale velocity and pressure fields
   - Set up ATD parameters

### Time-Stepping Loop

The main time-stepping loop performs the following operations for each time step:

1. **Save Previous Time Level**
   - Copy current solution (u, v, p, t) to previous time level (un, vn, pn, tn)
   - If small-scale enabled: save small-scale fields (usn, vsn, tsn)

2. **Momentum-Energy Iterations** (inner loop for coupling)
   - **Auxiliary Momentum Equations** (`nAuxMomentum` in `momentum.f`)
     - Newton-Kantorovich iterations (quasilinearization)
     - Solve X-momentum equation using Douglas-Gunn ADI
     - Solve Y-momentum equation using Douglas-Gunn ADI
     - Iterate until convergence or max iterations
   - **Shuman Filtering** (optional)
     - Apply filter to u-component if enabled
     - Apply filter to v-component if enabled
   - **Velocity Boundary Conditions**
     - Apply boundary conditions to auxiliary velocities
   - **Pressure Poisson Equation** (`Ppe` in `pressure.f`)
     - Compute divergence of auxiliary velocity
     - Solve using SOR (Successive Over-Relaxation)
     - Multiple solver options: point SOR, line SOR, red-black variants
   - **Projection Step** (`Project` in `utility.f`)
     - Project auxiliary velocity onto divergence-free space
     - Update velocities to time level n+1
   - **Thermal Energy Equation** (`ThermEnergy` in `thermal.f`, if enabled)
     - Solve temperature field using ADI method
   - **Equation of State** (`EqState` in `thermal.f`, if enabled)
     - Compute density from pressure and temperature
   - **Convergence Check**
     - Check difference between iterations
     - Exit if converged or max iterations reached

3. **Temperature Filtering** (if thermal enabled)
   - Apply Shuman filter to temperature field

4. **Small-Scale Computation** (if ATD enabled)
   - Subtract small-scale from previous time level
   - Compute new small-scale fluctuations (`SmallScale` in `small_scale.f`)
   - Add small-scale to large-scale fields

5. **Boundary Conditions**
   - Apply velocity boundary conditions
   - Apply pressure boundary conditions
   - Apply temperature boundary conditions (if enabled)

6. **Output Operations**
   - Time series output (if enabled)
   - Solution snapshots (if enabled)
   - Restart file writing (if enabled)
   - Particle trajectory tracking (if enabled)
   - Time-averaging accumulation (if enabled)

7. **Convergence/Divergence Check**
   - Check for solution divergence
   - Print difference norms if requested

## Key Algorithms

### Gresho's Projection Method

The projection method is used to enforce the incompressibility constraint. The algorithm consists of:

1. **Auxiliary Velocity Calculation**: Solve momentum equations without pressure gradient
2. **Pressure Poisson Equation**: Solve for pressure using divergence of auxiliary velocity
3. **Velocity Projection**: Correct velocities to be divergence-free

```fortran
! Pseudo-code
call nAuxMomentum(...)  ! Compute auxiliary velocities u*, v*
call Ppe(...)           ! Solve for pressure p
call Project(...)       ! Project velocities: u = u* - grad(p)
```

### Douglas-Gunn ADI Time-Splitting

The momentum equations are solved using ADI (Alternating Direction Implicit) method:

1. **X-Sweep**: Solve tridiagonal systems in i-direction
2. **Y-Sweep**: Solve tridiagonal systems in j-direction

This allows implicit treatment of diffusion terms while maintaining computational efficiency.

### Additive Turbulent Decomposition (ATD)

The ATD method decomposes the flow into large-scale and small-scale components:

- **Large-Scale**: Resolved by the numerical scheme
- **Small-Scale**: Modeled using linear combinations of chaotic maps

The small-scale components are computed separately and added to the large-scale solution.

### Successive Over-Relaxation (SOR)

The pressure Poisson equation is solved using iterative methods:

- **Point SOR**: Standard point-wise relaxation
- **Line SOR**: Line-wise relaxation (faster convergence)
- **Red-Black SOR**: Parallelizable variant
- **Red-Black Line SOR**: Combination of both approaches

## Data Structures

### Grid Variables

- `gx(0:mnx,0:mny)`, `gy(0:mnx,0:mny)` - Natural grid coordinates
- `ran, rbn, rgn` - Metric tensor components at natural points
- `rac, rbc, rgc` - Metric tensor components at cell centers
- `rau, rbu, rbv, rgv` - Metric tensor components at velocity points
- `dju, djv, djc, djn` - Jacobian determinants at different locations

### Solution Variables

- `u(0:mnx,0:mny)`, `v(0:mnx,0:mny)` - Velocity components
- `p(0:mnx,0:mny)` - Pressure
- `t(0:mnx,0:mny)` - Temperature
- `d(0:mnx,0:mny)` - Density
- `un, vn, pn, tn, dn` - Previous time level values
- `us, vs, ts` - Starred (intermediate) values

### Small-Scale Variables

- `uss(0:mnx,0:mny)`, `vss(0:mnx,0:mny)` - Small-scale velocities
- `pss(0:mnx,0:mny)` - Small-scale pressure
- `tss(0:mnx,0:mny)` - Small-scale temperature
- `usn, vsn, psn, tsn` - Previous time level small-scale values

### Boundary Condition Arrays

- `nReg(2)` - Number of regions in each direction
- `nRegBrd(mgri,mgrj,4)` - Region border locations
- `nRegType(mgri,mgrj)` - Region types (internal, porous, blockage)
- `nMomBdTp(mgri,mgrj,4)` - Momentum boundary types
- `nTemBdTp(mgri,mgrj,4)` - Thermal boundary types
- `dBCVal(mgri,mgrj,4,4)` - Boundary condition values

## Input/Output

### Input Files

- **input.dat**: Main input file with simulation parameters
- **bcond.dat**: Boundary condition file (defaults to input.dat)
- **grid file**: Grid coordinates (Plot3D format)
- **Restart file**: Previous solution (if restarting)

### Output Files

- **output.log**: Execution log file
- **output.xyz**: Final solution (Plot3D format)
- **output.qqq**: Additional solution fields
- **output.nam**: Variable names file
- **Snapshot files**: Time-dependent solution snapshots
- **Time series files**: Point-wise time histories
- **Restart files**: Checkpoint files for restarting
- **Trajectory files**: Particle position histories

## Module Descriptions

### main.f

Main program entry point. Contains:
- Variable declarations
- Input parameter reading
- Grid and boundary condition setup
- Time-stepping loop
- Output operations

**Key Subroutines Called:**
- `ReadParam` - Read input parameters
- `Grid` - Set up grid
- `SetUpBCs` - Set up boundary conditions
- `InitialCond` - Initialize solution
- `nAuxMomentum` - Solve momentum equations
- `Ppe` - Solve pressure equation
- `Project` - Project velocities
- `ThermEnergy` - Solve thermal equation
- `SmallScale` - Compute small-scale fields
- `SaveStdVarsP3D` - Save solution

### parse.f

Input parameter parsing module. Contains:
- `ReadParam` - Main parameter reading function
- `lReadIPFile` - Read input parameter file
- `lReadBCFile` - Read boundary condition file
- Token finding and parsing utilities

### grid.f

Grid processing module. Contains:
- `Grid` - Main grid setup driver
- `ReadGrid` - Read grid from file (Plot3D format)
- `MirrorPts` - Calculate mirror points
- `FullGrid` - Compute staggered grid coordinates
- `Metric` - Compute metric tensor components
- `CheckGridSize` - Validate grid dimensions

### bound_cond.f

Boundary condition module. Contains:
- `SetUpBCs` - Set up all boundary conditions
- `InitialCond` - Initialize or restart solution
- `VelBoundCond` - Apply velocity boundary conditions
- `PresBoundCond` - Apply pressure boundary conditions
- `TempBoundCond` - Apply temperature boundary conditions
- `SmlSclBC` - Apply small-scale boundary conditions
- `VelOutflowBCs` - Apply outflow boundary conditions

### momentum.f

Momentum equation solver. Contains:
- `nAuxMomentum` - Driver for auxiliary momentum equations
- `XMomentum` - X-momentum equation solver (ADI)
- `YMomentum` - Y-momentum equation solver (ADI)
- `ConvCoef` - Convective term discretization
- `DConvU`, `DConvV` - Convective term components
- `DDiffU`, `DDiffV` - Diffusive term components
- `PorosCoef` - Porous media terms
- `TridLU`, `AltTridLU` - Tridiagonal solvers

### pressure.f

Pressure Poisson equation solver. Contains:
- `Ppe` - Main pressure Poisson equation driver
- `Divergence` - Compute velocity divergence
- `RhsPpe` - Compute right-hand side of PPE
- `Sor` - Point SOR solver
- `SorRB` - Red-black point SOR
- `SorRBP` - Parallelizable red-black SOR
- `Slor` - Line SOR solver
- `SlorRB` - Red-black line SOR
- `SlorRBP` - Parallelizable red-black line SOR

### thermal.f

Thermal energy equation solver. Contains:
- `ThermEnergy` - Solve thermal energy equation
- `EqState` - Compute density from equation of state

### small_scale.f

ATD small-scale modeling. Contains:
- `SmallScale` - Compute small-scale fluctuations
- Uses chaotic maps to generate small-scale velocity components
- Enforces mass conservation for small-scale fields

### traject.f

Particle trajectory tracking. Contains:
- `InitTraject` - Initialize particle data
- `Traject` - Integrate particle positions
- `FwdEuler` - Forward Euler integration
- `HeunTrap` - Heun/Trapezoidal integration
- `TrajFunc` - Trajectory function evaluation
- `TrajJac` - Jacobian computation
- `FindPos` - Find grid indices for particle position
- `BiLinInterp` - Bilinear interpolation

### file_manip.f

File I/O operations. Contains:
- `SaveSolution` - Save solution in various formats
- `SaveGrid2DP3D` - Save grid coordinates
- `SaveRestart` - Save restart file
- `ReadRestart` - Read restart file
- `SaveFunc` - Save scalar function
- `SaveTraject` - Save particle trajectories
- `SaveTimeSrs` - Save time series data
- `SaveTmAvgP3D` - Save time-averaged quantities
- `SaveStdVarsP3D` - Save standard variables

### utility.f

Utility functions. Contains:
- `Filter` - Shuman filter for aliasing control
- `Project` - Velocity projection onto divergence-free space
- `PTDAvg` - Average pressure/temperature/density to natural points
- `VelAvg` - Average velocities to natural points
- `TAveraged` - Temperature averaging

### string.f

String manipulation utilities. Contains:
- `mToLowerC`, `mToUpperC` - Case conversion
- `mInitStr` - Initialize string
- `SyntaxError` - Error reporting
- `PrintTitle` - Print program title
- `PrintDiff` - Print difference norms
- `LogInputPar` - Log input parameters

## Configuration Parameters

Key parameters defined in `include/config.f`:

- `mnx`, `mny` - Maximum grid dimensions (default: 302)
- `mgri`, `mgrj` - Maximum number of regions (default: 20, 10)
- `mfnmlgth` - Maximum filename length (default: 32)
- `mlinelgt` - Maximum line length (default: 80)
- `mntr` - Maximum number of trajectories (default: 1)
- `maxtspts` - Maximum time-series points (default: 20)

## Compilation Directives

The program uses preprocessor directives for optional features:

- `_SMALLSCL_` - Enable ATD small-scale calculations
- `_TIMEAVG_` - Enable time-averaging mode
- `_TRAJECT_` - Enable particle trajectory tracking
- `_HP_` - HP-UX specific code (legacy)
- `DEBUG` - Additional debugging output

## See Also

- [Program Flow Diagram](diagrams/program-flow.md)
- [Data Flow Diagram](diagrams/data-flow.md)
- [Module Architecture Diagram](diagrams/module-architecture.md)
- [Time-Step Algorithm Diagram](diagrams/time-step-algorithm.md)
- [Porting Guide](gfortran-porting-guide.md)

