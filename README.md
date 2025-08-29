# `wolfd2`: 2D Incompressible Flow

`wolfd2` simulates 2D incompressible flow by numerically solving the
Navier-Stokes equations over generalized coordinates. Pressure-velocity coupling
is enforced using Gresho's 'Projection 1' method. Spatial derivatives are
discretized over finite volumes using cell-based indexing.  Time derivatives 
are integrated using Douglas-Gunn ADI time-splitting.  

The program also solves a thermal energy equation.  Buoyancy source terms are
included in the momentum equations. The modified Shuman filter is used to
control aliasing and cell-Reynolds (or cell-Peclet) number oscillations.

Turbulent flow is simulated using the Additive Turbulent Decomposition (ATD)
method of Hylin and McDonough.  This implementation of ATD uses linear
combinations of (appropriately-scaled) logistic maps to simulate the small-scale
velocity fluctuations.

The program simulates porous media flow by solving the momentum and thermal
energy equations with Darcy and Forschheimer terms.

Particle trajectories are calculated using a Lagrangian approach.

This program was written as part of masters thesis research under the
supervision of Prof. James M. McDonough, Computational Fluid Dynamics
Laboratory, Dept. of Mechanical Engineering, University of Kentucky.

**Version**: 0.3i  
**Original Platform**: HP-UX Fortran  
**Current Platform**: GNU Fortran (gfortran)  
**Last Updated**: August 2025

## Project Structure

The project has been modernized with a clean, organized structure:

```
wolfd2/
├── Makefile          # Top-level build system
├── src/              # Source files (*.f)
├── include/          # Header files (*.h, *.f)
├── build/            # Object files (*.o)
├── bin/              # Executables
├── docs/             # Documentation
└── scratch/          # Temporary files
```

## Building

### Prerequisites

- GNU Fortran (gfortran) compiler
- Make utility

### Quick Start

```bash
# Setup project directories
make setup

# Build executable
make

# Run the program
./bin/wolfd2
```

### Build Options

```bash
# Standard build (default)
make

# Debug build with full checking
make debug

# Optimized build for performance
make optimized

# Build with profiling support
make profile

# Test compilation of individual files
make test-compile

# Clean build artifacts
make clean

# Deep clean (remove all generated files)
make distclean
```

### Build Configurations

- **Standard**: `-O2 -Wall -std=legacy -cpp`
- **Debug**: `-g -O0 -Wall -std=legacy -cpp -fcheck=all`
- **Optimized**: `-O3 -march=native -ffast-math -std=legacy -cpp`
- **Profile**: `-O2 -pg -std=legacy -cpp`

## Porting Status

This codebase has been successfully ported from HP-UX Fortran to GNU Fortran (gfortran). Key changes include:

- ✅ HP-UX specific functions replaced with portable alternatives
- ✅ Compiler directives updated for gfortran compatibility
- ✅ Modern project structure implemented
- ✅ Comprehensive build system with multiple configurations

For detailed porting information, see `docs/gfortran-porting-guide.md`.

## Documentation

- `docs/gfortran-porting-guide.md` - Comprehensive porting guide and technical details
- `include/` - Header files and configuration
- `src/` - Source code with detailed comments

## License

This software was developed as part of academic research. Please contact the original authors for licensing information.
