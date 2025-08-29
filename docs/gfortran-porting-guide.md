# GNU Fortran Porting Guide for wolfd2

## Overview

This document outlines the process and requirements for porting the `wolfd2` Fortran 77 codebase from HP-UX Fortran to GNU Fortran (gfortran). The original code was written for HP-UX systems and uses several HP-UX-specific features that need to be replaced with portable alternatives.

## Project Background

`wolfd2` is a 2D incompressible flow simulation program that solves the Navier-Stokes equations using:
- Gresho's 'Projection 1' method for pressure-velocity coupling
- Finite volume discretization with cell-based indexing
- Douglas-Gunn ADI time-splitting for time integration
- Additive Turbulent Decomposition (ATD) for turbulence modeling
- Particle trajectory tracking capabilities

**Version**: 0.3i  
**Original Platform**: HP-UX Fortran 77  
**Target Platform**: GNU Fortran (gfortran)

## Key Issues Identified

### 1. HP-UX Specific Compiler Directives

The code uses several HP-UX specific compiler directives that must be replaced:

#### `$HP9000_800 INTRINSICS`
- **Location**: `src/grid.f` (line 145), `src/parse.h` (line 110)
- **Purpose**: Enables HP-UX specific intrinsic functions like `jnum()`, `dnum()`, and other HP-UX extensions
- **Action Required**: Remove these directives entirely
- **Impact**: This directive enables non-standard HP-UX intrinsic functions that are not available in gfortran

#### `$NOSTANDARD SYSTEM`
- **Location**: `src/main.f` (line 65)
- **Purpose**: Enables non-standard system functions like `secnds()` timing function
- **Action Required**: Remove this directive
- **Impact**: This directive allows the use of HP-UX specific system calls that are not portable

### 2. HP-UX Specific Functions

#### `secnds()` Timing Function
- **Location**: `src/main.f` (lines 305, 1356)
- **Usage**: 
  ```fortran
  rtsec = secnds(0.0)        ! Start timing
  rtsec = secnds(rtsec)      ! End timing and get elapsed time
  ```
- **Replacement**: Use standard Fortran `cpu_time()`
  ```fortran
  call cpu_time(rtsec_start)     ! Start timing
  call cpu_time(rtsec_end)       ! End timing
  rtsec = rtsec_end - rtsec_start ! Calculate elapsed time
  ```

#### `jnum()` and `dnum()` String Conversion Functions
- **Location**: `src/string.f` (lines 389, 415)
- **Usage**: Convert strings to integers and doubles
- **HP-UX Format**: Uses non-standard format specifiers `'(i)'` and `'(e)'`
- **Replacement**: Use standard Fortran `read` statements
  ```fortran
  ! Replace jnum(cNum) with:
  read(cNum,*) jnum
  
  ! Replace dnum(cNum) with:
  read(cNum,*) dnum
  ```

### 3. Conditional Compilation Blocks

The code extensively uses `#ifdef _HP_` conditional compilation:

#### Files with HP-UX Conditionals:
- `src/main.f` (multiple locations)
- `src/grid.f` (lines 143-147, 168)
- `src/string.f` (lines around 425)
- `src/parse.f` (lines 98, 1032, 2053)
- `src/parse.h` (lines 108-112)

#### Required Actions:
1. Remove `#ifdef _HP_` and `#endif _HP_` blocks
2. Keep the HP-UX specific code but make it portable
3. Update preprocessor definitions

#### Conditional Compilation Analysis:
The `#ifdef _HP_` blocks contain:
- HP-UX compiler directives (`$HP9000_800 INTRINSICS`, `$NOSTANDARD SYSTEM`)
- HP-UX specific variable declarations
- HP-UX intrinsic function calls
- Platform-specific code paths

### 4. Data Type Declarations

#### Current Usage:
```fortran
#define   REAL     real*8
#define   FLOAT    real*4
#define   INTEGER  integer*4
#define   LOGICAL  logical*4
```

#### Compatibility:
- **Status**: Generally compatible with gfortran
- **Recommendation**: Consider modern Fortran syntax for future maintenance
- **Modern Alternative**:
  ```fortran
  #define   REAL     real(kind=8)
  #define   FLOAT    real(kind=4)
  #define   INTEGER  integer(kind=4)
  #define   LOGICAL  logical(kind=4)
  ```

### 5. File I/O Operations

#### Status: Compatible
- Standard Fortran I/O operations used throughout
- No HP-UX specific I/O features identified
- Should work with gfortran without changes

## GNU Fortran Preprocessing and Compiler Directives

### Understanding gfortran Preprocessing

#### Preprocessing Capabilities
gfortran uses the C preprocessor (cpp) for handling preprocessor directives:

- **File Extensions**: `.F` and `.F90` files are automatically preprocessed
- **Manual Preprocessing**: Use `-cpp` flag to enable preprocessing for `.f` files
- **Preprocessor**: Uses the same C preprocessor as GCC (cpp)

#### Supported Preprocessor Directives
gfortran supports standard C preprocessor directives:
- `#include` - File inclusion
- `#define` - Macro definitions
- `#ifdef`, `#ifndef`, `#if`, `#else`, `#elif`, `#endif` - Conditional compilation
- `#undef` - Undefine macros
- `#error`, `#warning` - Error and warning directives

#### Preprocessing Commands
```bash
# Enable preprocessing for .f files
gfortran -cpp -c file.f

# Preprocess only (output to stdout)
gfortran -cpp -E file.f

# Preprocess with specific definitions
gfortran -cpp -D_HP_ -E file.f

# Pass options to preprocessor
gfortran -cpp -Wp,-I/path/to/include -c file.f
```

#### HP-UX vs gfortran Preprocessing Differences

| Feature | HP-UX Fortran | gfortran | Notes |
|---------|---------------|----------|-------|
| Preprocessor | Built-in | External (cpp) | gfortran uses C preprocessor |
| File Extensions | `.f` (with directives) | `.F` or `-cpp` flag | `.f` files need `-cpp` flag |
| Compiler Directives | `$DIRECTIVE` | Not supported | Must be removed |
| Conditional Compilation | `#ifdef` | `#ifdef` | Compatible |
| Include Files | `#include` | `#include` | Compatible |
| Macro Definitions | `#define` | `#define` | Compatible |

### HP-UX Compiler Directives Analysis

#### Non-Standard HP-UX Directives
HP-UX Fortran uses `$` prefixed directives that are not part of the Fortran standard:

1. **`$HP9000_800 INTRINSICS`**
   - **Purpose**: Enables HP-UX specific intrinsic functions
   - **Functions Enabled**: `jnum()`, `dnum()`, `secnds()`, and other HP-UX extensions
   - **gfortran Equivalent**: None - functions must be replaced with portable alternatives

2. **`$NOSTANDARD SYSTEM`**
   - **Purpose**: Enables non-standard system functions
   - **Functions Enabled**: `secnds()` timing function
   - **gfortran Equivalent**: Use `cpu_time()` intrinsic function

3. **`$HP9000_800`**
   - **Purpose**: Enables HP-UX 9000/800 series specific features
   - **gfortran Equivalent**: None - platform-specific optimizations

#### Standard Preprocessor Directives
These are compatible with gfortran:
- `#include "file.h"` - File inclusion
- `#define MACRO value` - Macro definitions
- `#ifdef _HP_` - Conditional compilation
- `#endif` - End conditional blocks

## Required Changes

#### Current HP-UX Makefile:
```makefile
F77 = f77
FFLAGS = +O3
FPPDEFS = -D_HP_ -D_SMALLSCL_ -D_TIMEAVG_
```

#### Recommended gfortran Makefile:
```makefile
F77 = gfortran
FFLAGS = -O3 -Wall -std=legacy
FPPDEFS = -D_SMALLSCL_ -D_TIMEAVG_
```

#### Key Changes:
- Replace `f77` with `gfortran`
- Replace HP-UX optimization flags (`+O3`) with gfortran equivalents (`-O3`)
- Add `-Wall` for warnings
- Add `-std=legacy` for Fortran 77 compatibility
- Remove `-D_HP_` from preprocessor definitions

#### Enhanced gfortran Makefile with Preprocessing:
```makefile
F77 = gfortran
FFLAGS = -O3 -Wall -std=legacy -cpp
FPPDEFS = -D_SMALLSCL_ -D_TIMEAVG_
INCLUDE = -I./src

# For debugging
FFLAGS_DEBUG = -g -O0 -Wall -std=legacy -cpp -fcheck=all
FPPDEFS_DEBUG = -D_SMALLSCL_ -D_TIMEAVG_ -DDEBUG

# For optimized compilation
FFLAGS_OPT = -O3 -march=native -ffast-math -std=legacy -cpp
FPPDEFS_OPT = -D_SMALLSCL_ -D_TIMEAVG_ -DNDEBUG
```

#### Preprocessing Considerations:
- **`-cpp` flag**: Enables preprocessing for `.f` files
- **`-I./src`**: Adds include path for header files
- **Conditional compilation**: Can use `-DDEBUG` for debug builds
- **File extensions**: Consider renaming `.f` to `.F` for automatic preprocessing

### 2. Replace HP-UX Functions

#### Timing Function Replacement:
```fortran
! In main.f, replace:
! rtsec = secnds(0.0)
call cpu_time(rtsec_start)

! And replace:
! rtsec = secnds(rtsec)
call cpu_time(rtsec_end)
rtsec = rtsec_end - rtsec_start
```

#### String Conversion Function Replacement:
```fortran
! In string.f, replace jnum function:
INTEGER function jnum(cNum)
  implicit none
  character*(*) cNum
  read(cNum,*) jnum
  return
end

! Replace dnum function:
REAL function dnum(cNum)
  implicit none
  character*(*) cNum
  read(cNum,*) dnum
  return
end
```

### 3. Remove HP-UX Compiler Directives

#### Remove from all files:
- `$HP9000_800 INTRINSICS`
- `$NOSTANDARD SYSTEM`

#### Update conditional compilation:
- Remove `#ifdef _HP_` and `#endif _HP_` blocks
- Keep the contained code but make it portable

#### Detailed Directive Migration Strategy:

##### Step 1: Identify and Remove HP-UX Directives
```bash
# Find all HP-UX directives in the codebase
grep -r "\$HP9000_800" src/
grep -r "\$NOSTANDARD" src/
```

##### Step 2: Handle Conditional Compilation Blocks
**Before (HP-UX):**
```fortran
#ifdef _HP_
$HP9000_800 INTRINSICS
#endif _HP_

#ifdef _HP_
$NOSTANDARD SYSTEM
#endif _HP_

#ifdef _HP_
      character*(mfnmlgth) wrkStr
#endif _HP_
```

**After (gfortran):**
```fortran
! Remove HP-UX directives entirely
! Keep portable variable declarations
      character*(mfnmlgth) wrkStr
```

##### Step 3: Update Preprocessor Definitions
**Remove from makefile:**
```makefile
# Remove this line:
FPPDEFS = -D_HP_ -D_SMALLSCL_ -D_TIMEAVG_

# Keep only portable definitions:
FPPDEFS = -D_SMALLSCL_ -D_TIMEAVG_
```

##### Step 4: Verify Preprocessing
```bash
# Test preprocessing without HP-UX definitions
gfortran -cpp -E src/main.f | grep -i "error\|warning"

# Test with specific definitions
gfortran -cpp -D_SMALLSCL_ -D_TIMEAVG_ -E src/main.f
```

### 4. Optional: Modernize Data Types

Consider updating type definitions for better portability:
```fortran
! In wolfd2.h, consider:
#define   REAL     real(kind=8)
#define   FLOAT    real(kind=4)
#define   INTEGER  integer(kind=4)
#define   LOGICAL  logical(kind=4)
```

## Additional Information Needed

### 1. Input/Output File Formats
- **Grid files**: Format and structure
- **Boundary condition files**: Parsing requirements
- **Parameter files**: Configuration format
- **Output files**: Expected formats

### 2. Dependencies
- Check for external library dependencies
- Verify mathematical library requirements
- Identify any system-specific calls

### 3. Test Cases
- Sample input files for testing
- Expected output files for validation
- Performance benchmarks

### 4. Performance Requirements
- Original HP-UX optimization levels
- Target performance expectations
- Memory usage requirements

## Porting Strategy

### Phase 1: Basic Compilation
1. Create gfortran makefile
2. Remove HP-UX compiler directives
3. Replace HP-UX specific functions
4. Test basic compilation

### Phase 2: Functionality Testing
1. Create test input files
2. Run basic functionality tests
3. Compare output with expected results
4. Debug any runtime issues

### Phase 3: Performance Optimization
1. Profile the application
2. Optimize gfortran compiler flags
3. Compare performance with original
4. Fine-tune as needed

### Phase 4: Validation
1. Run comprehensive test suite
2. Validate numerical accuracy
3. Performance benchmarking
4. Documentation updates

## Compiler Flags Recommendations

### Basic Compilation:
```bash
gfortran -O2 -Wall -std=legacy -cpp -o wolfd2 *.f
```

### Debug Compilation:
```bash
gfortran -g -O0 -Wall -std=legacy -cpp -fcheck=all -o wolfd2 *.f
```

### Optimized Compilation:
```bash
gfortran -O3 -march=native -ffast-math -std=legacy -cpp -o wolfd2 *.f
```

### Preprocessing-Specific Flags:
```bash
# Enable preprocessing for .f files
gfortran -cpp -c file.f

# Preprocess with specific definitions
gfortran -cpp -D_SMALLSCL_ -D_TIMEAVG_ -c file.f

# Preprocess only (for debugging)
gfortran -cpp -E file.f

# Pass options to preprocessor
gfortran -cpp -Wp,-I./src -c file.f
```

## Potential Issues and Solutions

### 1. Precision Differences
- **Issue**: Different floating-point precision between HP-UX and gfortran
- **Solution**: Use explicit kind specifications and test numerical accuracy

### 2. Endianness
- **Issue**: Different byte order for binary I/O files
- **Solution**: Test with sample files, convert if necessary

### 3. System Calls
- **Issue**: Any remaining system-specific calls
- **Solution**: Replace with portable alternatives

### 4. Performance
- **Issue**: Different optimization strategies
- **Solution**: Profile and tune gfortran flags

### 5. Preprocessing Issues
- **Issue**: HP-UX directives not recognized by gfortran
- **Solution**: Remove all `$` prefixed directives and replace with portable code
- **Issue**: Conditional compilation blocks may cause errors
- **Solution**: Remove `#ifdef _HP_` blocks and keep portable code

### 6. Compiler Directive Compatibility
- **Issue**: HP-UX specific intrinsics not available
- **Solution**: Replace with standard Fortran functions or custom implementations
- **Issue**: Non-standard format specifiers in I/O
- **Solution**: Use standard Fortran format specifiers

## Testing Checklist

- [ ] Compilation without errors
- [ ] No warnings (or acceptable warnings)
- [ ] Preprocessing works correctly
- [ ] All HP-UX directives removed
- [ ] Conditional compilation blocks handled
- [ ] Basic functionality test
- [ ] Input file parsing test
- [ ] Output file generation test
- [ ] Numerical accuracy validation
- [ ] Performance comparison
- [ ] Memory usage verification
- [ ] Preprocessor definitions working
- [ ] Include files processed correctly

## Conclusion

The port from HP-UX Fortran to GNU Fortran is feasible with the identified changes. The main challenges are:

1. **HP-UX specific functions** (`secnds`, `jnum`, `dnum`)
2. **Compiler directives** (`$HP9000_800 INTRINSICS`, `$NOSTANDARD SYSTEM`)
3. **Conditional compilation** blocks
4. **Preprocessing differences** between HP-UX and gfortran

The codebase uses standard Fortran 77 features for most operations, making it generally compatible with gfortran. The porting effort should focus on:

- **Replacing HP-UX specific features** with portable alternatives
- **Understanding gfortran preprocessing** capabilities and limitations
- **Properly handling compiler directives** and conditional compilation
- **Maintaining numerical accuracy** and performance
- **Testing preprocessing** and compilation thoroughly

### Key Success Factors:
- **Preprocessing**: Use `-cpp` flag or rename files to `.F` extension
- **Directive Migration**: Remove all HP-UX specific `$` directives
- **Function Replacement**: Replace HP-UX intrinsics with portable alternatives
- **Testing**: Comprehensive testing of preprocessing and compilation

## References

- [GNU Fortran Documentation](https://gcc.gnu.org/fortran/)
- [Fortran 77 Standard](https://wg5-fortran.org/f77-std.html)
- [Porting Guide for HP-UX to Linux](https://gcc.gnu.org/onlinedocs/gfortran/)

---

**Document Version**: 1.0  
**Last Updated**: [Current Date]  
**Status**: Draft for Review
