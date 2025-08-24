########################################################################
#  Generic description file for FORTRAN programs                       #
#  Version 1.2b                                                        #
#  vegar980814                                                         #
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
SHELL=/usr/bin/sh
#
.SUFFIXES: .f
#
########################################################################
# Objects of src files requiring preprocessing
#
OBJECTS = main.o grid.o momentum.o pressure.o thermal.o \
          smallScale.o traject.o miscel.o boundCond.o \
          fileManip.o string.o parse.o
#
HEADERS = config.f parse.h $(TARGET).h
SOURCE = $(OBJECTS:.o=.F)
FILES = makefile $(HEADERS) $(SOURCE) 
#
########################################################################
# Standard definitions
#
LIBDIRS =
LIBS =
INCLUDE =
F77 = f77
#
########################################################################
# Installation directory
#
INSTDIR = ../bin
#
########################################################################
# Compiler directives (HP)
#
#FFLAGS =
FFLAGS = +O3
#FFLAGS = +O1 +T 
#FFLAGS = +O1
#FFLAGS = +O2
#FFLAGS = +O2 +T
#FFLAGS = +O3 +T
#FFLAGS = -C
#FFLAGS = -C -g -y -L
#FFLAGS = -C -g
#FFLAGS = -C -g +T
#FFLAGS = -C -g +T -D
#FFLAGS = +O3 +ttyunbuf
#FFLAGS = +O3 +Onoparallel
#FFLAGS = +Oparallel +Oall
#
########################################################################
# Compiler directives for HP parallel servers:
#
#FFPAR1 = +O3 +U77 +Oparallel +Odataprefetch +Olibcalls +Onoautopar
#FFLAGS = $(FFPAR1) +Onofail_safe
#FFLAGS = $(FFPAR1) +pa +pal
#FFLAGS = $(FFPAR1) +Oreport=all +Oinfo +Onofail_safe
#LIBS = -Wl,+FPD -Wl,-a,archive_shared
#
########################################################################
# C preprocessor directive definitions
#
#FPPDEFS =
#FPPDEFS = -D_DATANH_
#FPPDEFS = -D_HP_
#FPPDEFS = -D_HP_ -D_SMALLSCL_
FPPDEFS = -D_HP_ -D_SMALLSCL_ -D_TIMEAVG_
#FPPDEFS = -D_HP_ -D_SMALLSCL_ -D_TRAJECT_ -D_TIMEAVG_
#
########################################################################
# Backup file declarations
#
DATESTR = $(HOME)/progs/tools/dateStr/dateStr -f
BAKDIR = ./bak/
#
BAKFILE = $(TARGET)-v$(VER)-`$(DATESTR)`.tar
# 
########################################################################
# `Update' file declarations
#
UPFLT = update-$(TARGET).tar
UPFLTG = $(UPFLT).gz
#
########################################################################
# Suffixes
#
.F.o :
	$(F77) $(FPPDEFS) $(FFLAGS) $(INCLUDE) -c $<

.f.o :
	$(F77) $(FFLAGS) $(INCLUDE) -c $<
#
########################################################################
# Main target
#
$(TARGET): $(OBJECTS)
	$(F77) $(FFLAGS) $(INCLUDE) $(LIBDIRS) $(LIBS) $(OBJECTS) -o $(TARGET)
	@size -v $(TARGET)
#
########################################################################
# Other targets
#
install: $(TARGET)
	@echo Copying executable to $(INSTDIR)
	@cp $(TARGET) $(INSTDIR)

########################################################################
backup:
	@echo Backup file time stamp: `$(DATESTR)`
	@tar cvf $(BAKFILE) $(FILES)
	@echo Compressing backup file...
	@gzip $(TARGET)-v$(VER)-*.tar
	@echo Moving compressed backup file to $(BAKDIR)
	@mv $(TARGET)-v$(VER)-*.tar.gz $(BAKDIR)

########################################################################
update:
	@if [ -f $(UPFLTG) ]; then \
	echo "Backing up current version of code."; make backup; \
	echo "Updating code from file $(UPFLTG)..."; \
	gunzip $(UPFLTG); tar xvf $(UPFLT); gzip $(UPFLT); \
	echo "Moving $(UPFLTG) to $(BAKDIR)..."; \
	mv $(UPFLTG) $(BKDIR); \
	else echo "Update file not found. Constructing one..."; \
	tar cvf $(UPFLT) $(SOURCE) $(HEADERS); gzip $(UPFLT); fi

########################################################################
cat:
	@echo 'Concatenating all relevant files into allFiles.txt'...
	@cat $(FILES) > allFiles.txt
	@echo Line count on allFiles.txt :  `wc -l allFiles.txt'

########################################################################
onefile:
	@echo 'Concatenating source files into '$(TARGET).F...
	@cat $(SOURCE) > $(TARGET).F
	@echo 'Line count on '$(TARGET).F :  `wc -l $(TARGET).F
	@echo 'Preprocessing '$(TARGET).F...
	@$(F77) -F $(TARGET).F
	@echo 'Line count on '$(TARGET).f :  `wc -l $(TARGET).f

########################################################################
clean:
	@echo 'Disk usage before: ' `du -s`
	@echo 'Cleaning...'
	@rm -f *.o $(TARGET)
	@rm -f allFiles.txt
	@rm -f a.out $(TARGET).f $(TARGET).F
	@echo 'Disk usage after : ' `du -s`

########################################################################
printout: $(SOURCE) $(HEADERS)
	cp $(SOURCE) $(HEADERS) ../doc/.
	(cd ../doc; src2tex $(SOURCE) $(HEADERS) > out.tex)
	(cd ../doc; latex out.tex; latex out.tex)
	(cd ../doc; dvips -o $(TARGET).ps out.dvi)

########################################################################
# Miscelaneous targets
#
fvpaths: fvpaths.f
	$(F77) $(FFLAGS) $(INCLUDE) -o fvpaths fvpaths.f

########################################################################
