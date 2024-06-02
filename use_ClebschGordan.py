from sympy import N
from sympy.physics.quantum.cg import CG, Wigner3j, Wigner6j

import ClebschGordan


def test_binomial():
    binomial = ClebschGordan.binomial_t()

    print(f"{binomial(5,2) = }")


def compare_sympy():
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

def calcCGs(N, _S, _Sprime, _Sdoubleprime):
    """Calculate Clebsch-Gordan coefficients for S x S' -> S''"

    Args:
        N (int): Specify Lie group SU(N)
        _S (list): N integers used to specify the irrep weight e.g. [2 1 0].
    """

    def init_weight(S, _S):
        "Assign weight element-wise. (Workaround for C++/swig)"
        for i, x in enumerate(_S):
            S[i+1] = x

    S = ClebschGordan.weight(N)
    Sprime = ClebschGordan.weight(N)
    Sdoubleprime = ClebschGordan.weight(N)

    #print('weight.dimension =', weight.dimension())

    init_weight(S, _S)
    init_weight(Sprime, _Sprime)
    init_weight(Sdoubleprime, _Sdoubleprime)
    dimS = S.dimension()
    dimSprime = Sprime.dimension()
    dimSdoubleprime = Sdoubleprime.dimension()

    coeffs = ClebschGordan.coefficients(Sdoubleprime, S, Sprime)
    for i in range(dimSdoubleprime):
        for j in range(dimSprime):
            for k in range(dimS):
                c = coeffs(j, k, 0, i)
                if abs(c) > 0.00001:
                    print(f"{j+1} {k+1} {i+1} {coeffs(j, k, 0, i)}")

if __name__ == "__main__":
    #calcCGs(2, [1, 0], [1, 0], [2, 0])
    calcCGs(3, [2, 1, 0], [2, 1, 0], [2, 1, 0])