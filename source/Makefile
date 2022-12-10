#-*- mode: makefile -*-
#-----------------------------------------------------------------------
#                 Intel Fortran Compiler (serial)
#-----------------------------------------------------------------------
# This file is part of the TEAKS package.
#
# Copyright (C) 2020-2026 Pai Li
# 2020-01-26  Pai Li
#-----------------------------------------------------------------------

#FC       = gfortran -c -g -pg
#LD       = gfortran -g -pg
FC       = ifort -c
LD       = ifort
DEBUG    =
FCFLAGS  = -O3 $(DEBUG)
LDFLAGS  = 
#LDFLAGS  = -static-libgfortran $(DEBUG)

OBJECTS  = io.o teio.o input.o species.o mesh.o \
           statistics.o gas.o core.o

TARGET   = teaks.x

#------------------------------- rules --------------------------------#

.SUFFIXES: .f90 .o .mod $(SUFFIXES)
.PHONY :  teaks  clean

teaks   : $(TARGET)

%.o : %.f90
	$(FC) $(FCFLAGS) $< -o $*.o

%.x : %.f90 $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $< $(OBJECTS)

clean :
	rm *.o *.mod teaks.x*

#----------------------- explicit dependencies ------------------------#

teio.o        : io.o
input.o       : io.o teio.o 
statistics.o  : io.o teio.o input.o
mesh.o        : io.o input.o statistics.o
species.o     : io.o teio.o input.o mesh.o statistics.o core.o
gas.o         : io.o teio.o input.o mesh.o species.o 
