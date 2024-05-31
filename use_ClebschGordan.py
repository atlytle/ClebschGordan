from sympy import N
from sympy.physics.quantum.cg import CG, Wigner3j, Wigner6j

import ClebschGordan

binomial = ClebschGordan.binomial_t()

#print(f"{binomial(5,2) = }")

#cg = CG(3/2, 3/2, 1/2, -1/2, 1, 1)
#print(f"{cg = }")
#print(f"{cg.doit() = }")
#print(f"{N(cg.doit()) = }")

# Let's compare some sympy results with ClebschGordan.

j1 = 1/2
j2 = 1/2

j3 = 1
print(f"{j3 = }")
for m1 in (-1/2, 1/2):
    for m2 in (-1/2, 1/2):
        for m3 in (1, 0, -1):
            cg = CG(j1, m1, j2, m2, j3, m3).doit()
            if cg != 0:
                print(f"{m1 = }, {m2 = }, {m3 = }")
                print(f"{cg = }")

print()
j3 = 0
print(f"{j3 = }")
for m1 in (-1/2, 1/2):
    for m2 in (-1/2, 1/2):
        for m3 in (0,):
            cg = CG(j1, m1, j2, m2, j3, m3).doit()
            if cg != 0:
                print(f"{m1 = }, {m2 = }, {m3 = }")
                print(f"{cg = }")
        
weight1 = ClebschGordan.weight(2)
weight2 = ClebschGordan.weight(2)
weight3 = ClebschGordan.weight(2)

#print('weight.dimension =', weight.dimension())
weight1[1] = 1
weight1[2] = 0
weight2[1] = 1
weight2[2] = 0
weight3[1] = 2
weight3[2] = 0
#print(f"{weight(1) = }")
coeffs = ClebschGordan.coefficients(weight3, weight1, weight2)
for i in range(2):
    for j in range(2):
        for k in range(2):
            print(f"{coeffs(i, j, 0, k)}")