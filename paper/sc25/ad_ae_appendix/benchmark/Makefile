FDPS_ROOT?=$(HOME)/FDPS
PS_PATH?=$(FDPS_ROOT)/src/
INC = -I$(PS_PATH)
CC = time g++
#CC = time mpicxx
#CFLAGS = -v
CFLAGS = -g -std=c++17 -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CFLAGS += -DPARTICLE_SIMULATOR_SPMOM_F32
CFLAGS += -DPARTICLE_SIMULATOR_USE_RADIX_SORT
#CFLAGS += -DPARTICLE_SIMULATOR_USE_SAMPLE_SORT
CFLAGS += -DPARTICLE_SIMULATOR_USE_64BIT_KEY
#CFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT
#use_phantom_grape_x86 = yes
#use_pikg_x86 = yes
#use_gpu_cuda = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

all:nbody.out

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(PS_PATH)/phantom_grape_x86/G5/newton/libpg5
INC += -I$(PG_ROOT)
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   = -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_BUILD =
PG_CLEAN = 
endif

ifeq ($(findstring yes, $(use_gpu_cuda) $(use_pikg_cuda)), yes)
#CUDA_ROOT ?= /usr/local/cuda
CUDA_ROOT ?= /usr/local/cuda-12.3
NVCC = time $(CUDA_ROOT)/bin/nvcc -O3 -std=c++17 -Xcompiler="$(CFLAGS)"
INC  += -I$(CUDA_ROOT)/samples/common/inc/
CFLAGS += -DENABLE_GPU_CUDA
CFLAGS += -DPIKG_MEASURE_CUDA_KERNEL_TIME
CLIBS = -L$(CUDA_ROOT)/lib64 -lcudart -lgomp

endif

ifeq ($(findstring yes, $(use_pikg_avx512) $(use_pikg_avx2) $(use_pikg_reference) $(use_pikg_cuda) $(use_pikg_mncore2)), yes)
PIKG_ROOT ?= ../../../../
PIKG = $(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference
INC += -I$(PIKG_ROOT)/inc
CFLAGS+= -DPIKG_USE_FDPS_VECTOR -DUSE_PIKG_KERNEL

ifeq ($(use_pikg_reference),yes)
endif

ifeq ($(use_pikg_avx512),yes)
CONVERSION_TYPE = AVX-512
CC=icpx
CFLAGS+= -xCORE_AVX512 -qopenmp
#CFLAGS+= -mavx512f -mavx512dq -ffast-math
endif

ifeq ($(use_pikg_avx2),yes)
CONVERSION_TYPE = AVX2
CC=icpx
CFLAGS+= -xCORE_AVX2 -qopenmp
endif

ifeq ($(use_pikg_cuda),yes)
CONVERSION_TYPE = CUDA
CFLAGS+=-DENABLE_PIKG_CUDA
endif

ifeq ($(use_pikg_mncore2),yes)
CONVERSION_TYPE = MNCore2
#PFCOMP_ROOT?=$(HOME)/work/pfcomp_latest
PFCOMP_ROOT?=/opt/pfn/pfcomp
MNCL_ROOT=$(PFCOMP_ROOT)/mncl
CC=g++
CFLAGS+= -DENABLE_PIKG_MNCORE2 -mavx2 -mfma -ffast-math #-march=x86-64-v4
CFLAGS+= -DONNX_ML=0
#CFLAGS+= -DPIKG_ENABLE_PERFETTO_TRACE
PKG_CONFIG_PATH=$(MNCL_ROOT)/lib/pkgconfig
CFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags mncl)
LDFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs mncl)
CFLAGS+=-I$(MNCL_ROOT)/include
LDFLAGS+=-lmncl -L$(MNCL_ROOT)/lib -Wl,-rpath=$(MNCL_ROOT)/lib
endif

PIKG_FLAGS= --conversion-type $(CONVERSION_TYPE)

HEDS+=kernel_pikg.hpp
kernel_pikg.hpp:	kernel_epep.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name EPIGrav --epj-name EPJGrav --force-name ForceGrav --kernel-name CalcGravityEpEp --class-file user-defined.hpp -i kernel_epep.pikg -o $@
else
PIKG_BUILD = 
endif

all:	nbody.out

ifeq ($(use_gpu_cuda), yes)
OBJS=force_gpu_cuda.o
force_gpu_cuda.o:	force_gpu_cuda.cu force_device.hpp
	$(NVCC) $(INC) -c -o $@ $<
endif

ifeq ($(use_pikg_cuda), yes)
OBJS=force_pikg_cuda.o
HEDS+=force_device.hpp
force_pikg_cuda.o:	force_pikg_device.cpp $(HEDS)
	$(NVCC) -gencode=arch=compute_90,code=compute_90 -x cu $(INC) -c -o $@ $<
endif
ifeq ($(use_pikg_mncore2),yes)
 OBJS=force_pikg_mncore2.o
HEDS+=force_device.hpp
force_pikg_mncore2.o:	force_pikg_device.cpp $(HEDS)
	$(CC) $(INC) $(CFLAGS) -c -o $@ $< $(CLIBS) $(LDFLAGS)
endif

nbody.out:	nbody.cpp $(OBJS) $(HEDS)
	$(PG_BUILD)
	$(CC) $(INC) $(CFLAGS) -o $@ $< $(OBJS) $(CLIBS) $(LDFLAGS)

clean:
	rm -f *.o nbody.out kernel_pikg.hpp *~

distclean: clean
	$(PG_CLEAN)
	rm -f nbody.out
	rm -rf result

test: 
	# This command is only for FDPS developers.
	./test.py

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
