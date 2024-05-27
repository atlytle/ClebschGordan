all: ClebschGordan

swig: ClebschGordan.py ClebschGordan_wrap.cxx

ClebschGordan: ClebschGordan.cpp
	g++ ClebschGordan.cpp -o ClebschGordan -llapack

ClebschGordan.py ClebschGordan_wrap.cxx: ClebschGordan.i ClebschGordan.h ClebschGordan.c
	swig -python ClebschGordan.i
	python setup_ClebschGordan.py build_ext --inplace

clean:
	rm ClebschGordan
	rm ClebschGordan.py ClebschGordan_wrap.cxx _Clebsch*
	rm -rf build/
