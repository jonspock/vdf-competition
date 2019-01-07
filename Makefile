
build = build

tb: $(build)/Makefile
	$(MAKE) -C $(build)

$(build)/Makefile: CMakeLists.txt
	mkdir -p $(build) && cd $(build) && cmake ..

update_ignore:
	file -F" " -i * | grep exec | awk '{print $$1}' | xargs -0 >.gitignore

clean:
	# remove executables in current directory
	file -F" " -i * | grep exec | awk '{print $$1}' | xargs rm -f
	rm -rf *.a *.so $(build) test_*
