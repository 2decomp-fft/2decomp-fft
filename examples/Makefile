.PHONY: test2d timing halo_test 

# Just build the examples
examples: test2d timing halo_test 
	@echo "Built the examples"

test2d:
	$(MAKE) -C $@ $@
timing:
	$(MAKE) -C $@ $@
halo_test:
	$(MAKE) -C $@ $@

check:
	cd test2d; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@

export
