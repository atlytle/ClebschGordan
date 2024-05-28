cpp: ClebschGordan

swig: ClebschGordan.py _ClebschGordan.so

ClebschGordan: ClebschGordan.cpp
	g++ ClebschGordan.cpp -o ClebschGordan -llapack

ClebschGordan.py ClebschGordan_wrap.cxx: ClebschGordan.i ClebschGordan.h ClebschGordan.cpp
	swig -c++ -python ClebschGordan.i

ClebschGordan.o ClebschGordan_wrap.o: ClebschGordan.cpp ClebschGordan_wrap.cxx	
	g++ -fPIC -c ClebschGordan.cpp ClebschGordan_wrap.cxx -I/usr/include/python3.11 -llapack

_ClebschGordan.so: ClebschGordan.o ClebschGordan_wrap.o
	g++ -shared ClebschGordan.o ClebschGordan_wrap.o -llapack -o _ClebschGordan.so 

clean:
	rm -f ClebschGordan
	rm -f ClebschGordan.py ClebschGordan_wrap.cxx _Clebsch* *.o
	#rm -rf build/
