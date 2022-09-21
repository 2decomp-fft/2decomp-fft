.PHONY: test2d fft_test_c2c fft_test_r2c timing halo_test io_test

# Just build the examples
examples: test2d fft_test_c2c fft_test_r2c timing halo_test io_test
	@echo "Built the examples"

test2d:
	$(MAKE) -C $@ $@
fft_test_c2c:
	$(MAKE) -C $@ $@
fft_test_r2c:
	$(MAKE) -C $@ $@
timing:
	$(MAKE) -C $@ $@
halo_test:
	$(MAKE) -C $@ $@
io_test:
	$(MAKE) -C $@ all

check:
	cd test2d; $(MAKE) $@
	cd fft_test_c2c; $(MAKE) $@
	cd fft_test_r2c; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

check:
	cd test2d; $(MAKE) $@
	cd fft_test_c2c; $(MAKE) $@
	#cd fft_test_r2c; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd fft_test_c2c; $(MAKE) $@
	cd fft_test_r2c; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@

export
