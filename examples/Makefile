.PHONY: test2d fft_physical_z fft_physical_x halo_test init_test io_test 

# Just build the examples
examples: test2d fft_physical_z fft_physical_x halo_test init_test io_test
	@echo "Built the examples"

test2d:
	$(MAKE) -C $@ $@
fft_physical_z:
	$(MAKE) -C $@ $@
fft_physical_x:
	$(MAKE) -C $@ all
halo_test:
	$(MAKE) -C $@ $@
init_test:
	$(MAKE) -C $@ all
io_test:
	$(MAKE) -C $@ all

check:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd init_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd init_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

export
