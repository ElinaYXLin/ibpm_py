# build/Makefile
#
# Compiles src/*.cc into libibpm.a, and links the driver programs
# (ibpm, checkgeom) against it + FFTW3.
#
# NOTE: this file was reconstructed -- it was missing from this checkout
# (build/ is gitignored and, unusually, the Makefile that's normally
# committed there was absent too). Layout follows test/Makefile's
# expectations (BUILDDIR=../build, target libibpm.a) and src/'s own
# file layout (two driver .cc files with main(), everything else is
# library code).

include ../config/make.inc

SRCDIR = ../src
CPPFLAGS += $(include_dirs) -I$(SRCDIR)
LDFLAGS += $(lib_dirs)
LIBS = -lfftw3 -lm

IBPMLIB = libibpm.a

DRIVERS = ibpm.cc checkgeom.cc
ALL_SRCS = $(notdir $(wildcard $(SRCDIR)/*.cc))
LIB_SRCS = $(filter-out $(DRIVERS), $(ALL_SRCS))
LIB_OBJS = $(LIB_SRCS:.cc=.o)

.PHONY: all clean distclean

all: ibpm checkgeom

$(IBPMLIB): $(LIB_OBJS)
	ar rcs $@ $(LIB_OBJS)

ibpm: ibpm.o $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ ibpm.o $(IBPMLIB) $(LIBS)

checkgeom: checkgeom.o $(IBPMLIB)
	$(CXX) $(LDFLAGS) -o $@ checkgeom.o $(IBPMLIB) $(LIBS)

%.o : $(SRCDIR)/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o $(IBPMLIB)

distclean: clean
	rm -f ibpm checkgeom
