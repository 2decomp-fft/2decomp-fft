.PHONY: init_test test2d fft_physical_x fft_physical_z halo_test io_test 

# Just build the examples
examples: init_test test2d fft_physical_x fft_physical_z halo_test io_test
	@echo "Built the examples"

init_test:
	$(MAKE) -C $@ all
test2d:
	$(MAKE) -C $@ all
fft_physical_x:
	$(MAKE) -C $@ all
fft_physical_z:
	$(MAKE) -C $@ all
halo_test:
	$(MAKE) -C $@ $@
io_test:
	$(MAKE) -C $@ all

check:
	cd init_test; $(MAKE) $@
	cd test2d; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

clean:
	cd init_test; $(MAKE) $@
	cd test2d; $(MAKE) $@
	cd fft_physical_x; $(MAKE) $@
	cd fft_physical_z; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

export
