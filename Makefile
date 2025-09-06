########################################################################
#  GNU Fortran Makefile for wolfd2                                       #
#  Modern Project Structure: build/, bin/, include/, src/               #
#  Ported from HP-UX Fortran to GNU Fortran (gfortran)                  #
########################################################################
#                                                                      #
########################################################################
# Target name and revision number
#
TARGET = wolfd2
VER = 0.3i
#
########################################################################
#
SHELL=/bin/sh
#
.SUFFIXES: .f
#
########################################################################
# Directory structure
#
SRCDIR = src
BUILDDIR = build
BINDIR = bin
INCLUDEDIR = include
#
########################################################################
# Objects of src files requiring preprocessing
#
SOURCES = $(wildcard $(SRCDIR)/*.f)
OBJECTS = $(SOURCES:$(SRCDIR)/%.f=$(BUILDDIR)/%.o)
HEADERS = $(wildcard $(INCLUDEDIR)/*.h $(INCLUDEDIR)/*.f)
#
########################################################################
# Standard definitions
#
LIBDIRS =
LIBS =
INCLUDE = -I$(INCLUDEDIR)
F77 = gfortran
#
########################################################################
# Compiler directives (GNU Fortran)
#
# Basic compilation
FFLAGS = -O2 -Wall -std=legacy -cpp
FPPDEFS = -D_SMALLSCL_ -D_TIMEAVG_ -D_TRAJECT_

# Debug compilation
FFLAGS_DEBUG = -g -O0 -Wall -std=legacy -cpp -fcheck=all
FPPDEFS_DEBUG = -D_SMALLSCL_ -D_TIMEAVG_ -D_TRAJECT_ -DDEBUG

# Optimized compilation
FFLAGS_OPT = -O3 -march=native -ffast-math -std=legacy -cpp
FPPDEFS_OPT = -D_SMALLSCL_ -D_TIMEAVG_ -D_TRAJECT_ -DNDEBUG

# Profile compilation
FFLAGS_PROF = -O2 -pg -std=legacy -cpp
FPPDEFS_PROF = -D_SMALLSCL_ -D_TIMEAVG_ -D_TRAJECT_
#
########################################################################
# Suffixes
#
$(BUILDDIR)/%.o: $(SRCDIR)/%.f
	@mkdir -p $(BUILDDIR)
	$(F77) $(FPPDEFS) $(FFLAGS) $(INCLUDE) -c $< -o $@

$(BUILDDIR)/%_debug.o: $(SRCDIR)/%.f
	@mkdir -p $(BUILDDIR)
	$(F77) $(FPPDEFS_DEBUG) $(FFLAGS_DEBUG) $(INCLUDE) -c $< -o $@

$(BUILDDIR)/%_opt.o: $(SRCDIR)/%.f
	@mkdir -p $(BUILDDIR)
	$(F77) $(FPPDEFS_OPT) $(FFLAGS_OPT) $(INCLUDE) -c $< -o $@

$(BUILDDIR)/%_prof.o: $(SRCDIR)/%.f
	@mkdir -p $(BUILDDIR)
	$(F77) $(FPPDEFS_PROF) $(FFLAGS_PROF) $(INCLUDE) -c $< -o $@
#
########################################################################
# Main targets
#
.PHONY: all clean distclean test-compile preprocess-test help install debug optimized profile

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(F77) $(FFLAGS) $(INCLUDE) $(LIBDIRS) $(LIBS) $(OBJECTS) -o $@
	@echo "Build completed successfully"
	@ls -la $@

debug: $(BINDIR)/$(TARGET)_debug

$(BINDIR)/$(TARGET)_debug: $(OBJECTS:%.o=%_debug.o)
	@mkdir -p $(BINDIR)
	$(F77) $(FFLAGS_DEBUG) $(INCLUDE) $(LIBDIRS) $(LIBS) $^ -o $@
	@echo "Debug build completed successfully"

optimized: $(BINDIR)/$(TARGET)_opt

$(BINDIR)/$(TARGET)_opt: $(OBJECTS:%.o=%_opt.o)
	@mkdir -p $(BINDIR)
	$(F77) $(FFLAGS_OPT) $(INCLUDE) $(LIBDIRS) $(LIBS) $^ -o $@
	@echo "Optimized build completed successfully"

profile: $(BINDIR)/$(TARGET)_prof

$(BINDIR)/$(TARGET)_prof: $(OBJECTS:%.o=%_prof.o)
	@mkdir -p $(BINDIR)
	$(F77) $(FFLAGS_PROF) $(INCLUDE) $(LIBDIRS) $(LIBS) $^ -o $@
	@echo "Profile build completed successfully"
#
########################################################################
# Other targets
#
install: $(BINDIR)/$(TARGET)
	@echo "Executable is ready in $(BINDIR)/$(TARGET)"

########################################################################
clean:
	@echo 'Cleaning build artifacts...'
	@rm -rf $(BUILDDIR)/*.o $(BINDIR)/$(TARGET)*
	@rm -f $(BUILDDIR)/*.mod
	@echo 'Clean completed'

########################################################################
distclean: clean
	@echo 'Deep cleaning...'
	@rm -rf $(BUILDDIR) $(BINDIR)
	@rm -f allFiles.txt
	@rm -f a.out $(TARGET).f $(TARGET).F
	@rm -f gmon.out
	@echo 'Deep clean completed'

########################################################################
test-compile:
	@echo 'Testing compilation of individual files...'
	@mkdir -p $(BUILDDIR)
	@for file in $(SOURCES); do \
		echo "Compiling $$file..."; \
		$(F77) $(FPPDEFS) $(FFLAGS) $(INCLUDE) -c $$file -o $(BUILDDIR)/$$(basename $$file .f).o || exit 1; \
	done
	@echo 'Individual file compilation test completed'

########################################################################
preprocess-test:
	@echo 'Testing preprocessing...'
	@for file in $(SOURCES); do \
		echo "Preprocessing $$file..."; \
		$(F77) $(FPPDEFS) $(FFLAGS) $(INCLUDE) -E $$file > /dev/null || exit 1; \
	done
	@echo 'Preprocessing test completed'

########################################################################
setup:
	@echo 'Setting up project directories...'
	@mkdir -p $(BUILDDIR) $(BINDIR) $(INCLUDEDIR)
	@echo 'Project structure created:'
	@echo '  $(SRCDIR)/     - Source files'
	@echo '  $(BUILDDIR)/   - Object files'
	@echo '  $(BINDIR)/     - Executables'
	@echo '  $(INCLUDEDIR)/ - Header files'

########################################################################
help:
	@echo 'Available targets:'
	@echo '  all         - Build with standard optimization (default)'
	@echo '  debug       - Build with debug flags'
	@echo '  optimized   - Build with maximum optimization'
	@echo '  profile     - Build with profiling support'
	@echo '  install     - Install executable to $(BINDIR)'
	@echo '  clean       - Remove object files and executables'
	@echo '  distclean   - Deep clean (remove all generated files)'
	@echo '  test-compile - Test compilation of individual files'
	@echo '  preprocess-test - Test preprocessing of all files'
	@echo '  setup       - Create project directory structure'
	@echo '  help        - Show this help message'
	@echo ''
	@echo 'Project structure:'
	@echo '  $(SRCDIR)/     - Source files (*.f)'
	@echo '  $(BUILDDIR)/   - Object files (*.o)'
	@echo '  $(BINDIR)/     - Executables'
	@echo '  $(INCLUDEDIR)/ - Header files (*.h, *.f)'
	@echo ''
	@echo 'Compiler flags:'
	@echo '  Standard: $(FFLAGS)'
	@echo '  Debug: $(FFLAGS_DEBUG)'
	@echo '  Optimized: $(FFLAGS_OPT)'
	@echo '  Profile: $(FFLAGS_PROF)'

########################################################################
