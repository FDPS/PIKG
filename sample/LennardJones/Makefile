CC=g++
CCFLAGS= -O2 

RUBY=/usr/local/opt/ruby/bin/ruby
PIKG_ROOT=../../
PIKG=$(RUBY) -I$(PIKG_ROOT)/src $(PIKG_ROOT)/src/parserdriver.rb
CONVERSION_TYPE=reference
CCFLAGS += -I../../inc

PIKG_OPTION= --conversion-type $(CONVERSION_TYPE) --epi-name EP --epj-name EP --force-name Force
# AVX2 mode
#CONVERSION_TYPE=AVX2 
#CCFLAGS += -mavx2 -mfma -ffast-math

# AVX-512 mode
#CONVERSION_TYPE=AVX-512
#CCFLAGS += -mavx512f -mavx512dq -ffast-math

# ARM SVE mode
#CONVERSION_TYPE=A64FX
#CCFLAGS += -Kfast

all:	main.out
kernel.hpp:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
main.out:	main.cpp kernel.hpp Makefile
	$(CC) $(CCFLAGS) main.cpp -o $@
test:	main.out
	./main.out
clean:
		rm kernel.hpp main.out
