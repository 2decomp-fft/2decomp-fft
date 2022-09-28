.PHONY: test2d fft_physical_z fft_physical_x halo_test 

# Just build the examples
examples: test2d fft_physical_z fft_physical_x halo_test 
	@echo "Built the examples"

test2d:
	$(MAKE) -C $@ $@
fft_physical_z:
	$(MAKE) -C $@ $@
fft_physical_x:
	$(MAKE) -C $@ $@
halo_test:
	$(MAKE) -C $@ $@

check:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

export
