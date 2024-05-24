all: ClebschGordan

swig: example.py example_wrap.c

ClebschGordan: ClebschGordan.cpp
	g++ ClebschGordan.cpp -o ClebschGordan -llapack

example.py example_wrap.c: example.i example.h example.c
	swig -python example.i
	python setup_swig.py build_ext --inplace

clean:
	rm example.py example_wrap.c ClebschGordan
