.PHONY: test2d fft_physical_z halo_test 

# Just build the examples
examples: test2d fft_physical_z halo_test 
	@echo "Built the examples"

test2d:
	$(MAKE) -C $@ $@
fft_physical_z:
	$(MAKE) -C $@ $@
halo_test:
	$(MAKE) -C $@ $@

check:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd halo_test; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd halo_test; $(MAKE) $@

export
