import math

#from sympy import N
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
        _S (list): N integers to specify the irrep weight e.g. [2 1 0].
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
    print("Q' Q Q''")
    for i in range(dimSdoubleprime):
        for j in range(dimSprime):
            for k in range(dimS):
                c = coeffs(j, k, 0, i)
                if abs(c) > 0.00001:
                    print(f"{j+1} {k+1} {i +1} {coeffs(j, k, 0, i)}")

def CG_swig(N, w1, m1, w2, m2, w3, m3):
    """CG coefficient appearing in w1 x w2-> w3.

    Args:
        N (int): Specify Lie group SU(N)
        w1,... (list): N integers to specify irrep weight e.g. [2 1 0]
        m1,... (int): Quantum number associated to w1
    """

    def init_weight(S, _S):
        "Assign weight element-wise. (Workaround for C++/swig)"
        for i, x in enumerate(_S):
            S[i+1] = x

    S = ClebschGordan.weight(N)
    Sprime = ClebschGordan.weight(N)
    Sdoubleprime = ClebschGordan.weight(N)

    #print('weight.dimension =', weight.dimension())

    init_weight(S, w1)
    init_weight(Sprime, w2)
    init_weight(Sdoubleprime, w3)
    dimS = S.dimension()
    dimSprime = Sprime.dimension()
    dimSdoubleprime = Sdoubleprime.dimension()

    # Calculates all coefficients.
    coeffs = ClebschGordan.coefficients(Sdoubleprime, S, Sprime)

    return coeffs(m2, m1, 0, m3)  # Set multiplicity index to 0 for now..

def test_CG_swig():
    """Compare output of CG_swig with sympy.CG (SU(2) only).
    """
    N = 2
    j1 = 1/2
    j2 = 1/2
    j3 = 1
    w1 = [1, 0]
    w2 = [1, 0]
    w3 = [2, 0]

    m1 = -1/2
    m2 = -1/2
    m3 = -1

    cg = CG(j1, m1, j2, m2, j3, m3).doit()
    cg2 = CG_swig(N, w1, 0, w2, 0, w3, 0)
    print(cg)
    print(cg2)
    print(math.isclose(cg, cg2))


if __name__ == "__main__":
    #calcCGs(2, [1, 0], [1, 0], [2, 0])
    
    #calcCGs(3, [2, 1, 0], [2, 1, 0], [2, 1, 0])
    #print()
    #print(CG_swig(3, [2, 1, 0], 2, [2, 1, 0], 0, [2, 1, 0], 0))
    #print(CG_swig(3, [2, 1, 0], 0, [2, 1, 0], 2, [2, 1, 0], 0))
    test_CG_swig()