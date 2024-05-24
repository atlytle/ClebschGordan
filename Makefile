all: ClebschGordan

ClebschGordan: ClebschGordan.cpp
	g++ ClebschGordan.cpp -o ClebschGordan -llapack
