# `wolfd2` 2D Incompressible Flow

## Project Background

`wolfd2` is a 2D incompressible flow simulation program that solves the Navier-Stokes equations using:
- Gresho's 'Projection 1' method for pressure-velocity coupling
- Finite volume discretization with cell-based indexing
- Douglas-Gunn ADI time-splitting for time integration
- Additive Turbulent Decomposition (ATD) for turbulence modeling
- Lagrangian particle trajectory tracking

The program was written as part of master's thesis research under the
supervision of Prof. James M. McDonough, Computational Fluid Dynamics
Laboratory, Dept. of Mechanical Engineering, University of Kentucky.

**Version**: 0.3i  
**Original Platform**: HP-UX Fortran  
**Current Platform**: GNU Fortran (gfortran)

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
└── examples/         # Simulation examples
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

### Preprocessing directives

- **`D_SMALLSCL_`**: ATD small scale subroutines
- **`-D_TIMEAVG_`**: Time averaging subroutines
- **`-D_TRAJECT_`**: Lagrangian trajectory tracking
- **`-DDEBUG`**: Additional debugging output

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
