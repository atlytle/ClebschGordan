all: example.py example_wrap.c

example.py example_wrap.c: example.i example.h example.c
	swig -python example.i
	python setup_swig.py build_ext --inplace

clean:
	rm example.py example_wrap.c _example*
	rm -rf build/
