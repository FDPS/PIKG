PIKG_ROOT=../../../../
PIKG=$(PIKG_ROOT)/bin/pikg
PIKG_FLAGS= --epi-name EPI --epj-name EPJ --force-name FORCE
CONVERSION_TYPE?=reference
all:	gravity_reference.hpp gravity_avx2.hpp lennard_jones_reference.hpp lennard_jones_avx2.hpp

gravity_reference.hpp:	gravity.pikg class_gravity.hpp
	$(PIKG) $(PIKG_FLAGS) --class-file class_gravity.hpp --conversion-type reference -i $< -o $@
lennard_jones_reference.hpp:	lennard_jones.pikg class_lennard_jones.hpp
	$(PIKG) $(PIKG_FLAGS) --class-file class_lennard_jones.hpp --conversion-type reference -i $< -o $@
gravity_avx2.hpp:	gravity.pikg class_gravity.hpp
	$(PIKG) $(PIKG_FLAGS) --class-file class_gravity.hpp --conversion-type AVX2 -i $< -o $@
lennard_jones_avx2.hpp:	lennard_jones.pikg class_lennard_jones.hpp
	$(PIKG) $(PIKG_FLAGS) --class-file class_lennard_jones.hpp --conversion-type AVX2 -i $< -o $@
clean:
	rm -rf gravity*.hpp lennard_jones*.hpp
