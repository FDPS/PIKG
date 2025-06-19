PIKG_ROOT=../..
PIKG=$(PIKG_ROOT)/bin/pikg

CONVERSION_TYPE?=reference

CC=g++
CCFLAGS= -std=c++20 -O2 -g
INC=-I$(PIKG_ROOT)/inc
ifeq ($(CONVERSION_TYPE),CUDA)
  CUDA_ROOT=/usr/local/cuda
  CC=$(CUDA_ROOT)/bin/nvcc
  CCFLAGS= -std=c++20 -x cu -Xcompiler="-std=c++20 -O3 -DMULTI_WALK"
  INC+=-I/usr/local/cuda/include
endif
ifeq ($(CONVERSION_TYPE),AVX2)
  CCFLAGS+= -mavx2
endif
ifeq ($(CONVERSION_TYPE),AVX-512)
  CCFLAGS+= -mavx512f
endif

PIKGFLAGS=--conversion-type $(CONVERSION_TYPE) --kernel-name CalcHydroForceEpEp --class-file user_defined.h --epi-name EP_hydro --epj-name EP_hydro --force-name Force_hydro 

ifeq ($(CONVERSION_TYPE), A64FX)
  PIKGFLAGS+= --strip-mining 64
endif

all:	main

main:	main.cpp hydro_force.hpp
	$(CC) $(CCFLAGS) $(INC) $< -o $@

hydro_force.hpp:	hydro_force.pikg user_defined.h
	$(PIKG) -o $@ -i $< $(PIKGFLAGS)

clean:
	rm -f hydro_force.hpp main
