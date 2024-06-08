#!/usr/bin/env python3
# Andrew Lytle
# June 2024

import math
from itertools import product

import numpy as np

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
    # dimS = S.dimension()
    # dimSprime = Sprime.dimension()
    # dimSdoubleprime = Sdoubleprime.dimension()

    # Calculates all coefficients.
    coeffs = ClebschGordan.coefficients(Sdoubleprime, S, Sprime)

    return coeffs(m1, m2, 0, m3)  # Set multiplicity index to 0 for now..

def test_CG_swig(w1, w2, w3):
    """Compare output of CG_swig with sympy.CG (SU(2) only).
    """
    N = 2
    #w1 = [1, 0]
    #w2 = [1, 0]
    #w3 = [2, 0]
    j1 = w1[0]/2
    j2 = w2[0]/2
    j3 = w3[0]/2

    result = []
    for a1, a2, a3 in product(range(w1[0]+1), 
                              range(w2[0]+1), 
                              range(w3[0]+1)):
        m1 = -j1 + a1
        m2 = -j2 + a2
        m3 = -j3 + a3
        #print(f"{m1 = } {m2 = } {m3 = }")
        cg = CG(j1, m1, j2, m2, j3, m3).doit()
        cg2 = CG_swig(N, w1, a1, w2, a2, w3, a3)
        #print(cg)
        #print(cg2)
        _result = math.isclose(cg, cg2)
        #print(_result)
        result.append(_result)
    return np.array(result).all()

def compute_9R(A, B, C, D, E):
    "Implement Eq. 9 of [2101.10227]."

    def init_weight(S, _S):
        "Assign weight element-wise. (Workaround for C++/swig)"
        for i, x in enumerate(_S):
            S[i+1] = x

    result = 0
    N = 3  # SU(3).
    three = [1, 0, 0]  # Weight associated to the fundamental(3) of SU(3)
    _three = ClebschGordan.weight(N) # Fix, this can be 3bar too..
    init_weight(_three, three)
    dthree = _three.dimension()
    # Need to fix this sort of pattern..
    _A = ClebschGordan.weight(N)
    init_weight(_A, A)
    dA = _A.dimension()
    _B = ClebschGordan.weight(N)
    init_weight(_B, B)
    dB = _B.dimension()
    _C = ClebschGordan.weight(N)
    init_weight(_C, C)
    dC = _C.dimension()
    _D = ClebschGordan.weight(N)
    init_weight(_D, D)
    dD = _D.dimension()
    _E = ClebschGordan.weight(N)
    init_weight(_E, E)
    dE = _E.dimension()
    print(f"{dthree = } {dA = } {dB = } {dC = } {dD = } {dE = } ")
    result = 0
    for a, b, c, d, e, m in product(range(dA), range(dB), range(dC),
                                    range(dD), range(dE), range(dthree)):
        print(a,b,c,d,e,m)
        result += CG_swig(N, D, d, B, b, E, e)*CG_swig(N, A, a, B, b, C, c)*\
                  CG_swig(N, A, a, three, m, D, d)*CG_swig(N, C, c, three, m, E, e)
        print(CG_swig(N, D, d, B, b, E, e))
        print(CG_swig(N, D, d, B, b, E, e))
        print(CG_swig(N, A, a, three, m, D, d))
        print(CG_swig(N, C, c, three, m, E, e))
    return result

def debug_CG():
    N = 3
    print(CG_swig(N, [0, 0, 0], 0, [1, 0, 0], 1, [1, 0, 0], 0))

if __name__ == "__main__":
    #calcCGs(2, [1, 0], [1, 0], [2, 0])
    
    #calcCGs(3, [2, 1, 0], [2, 1, 0], [2, 1, 0])
    #print()
    #print(CG_swig(3, [2, 1, 0], 2, [2, 1, 0], 0, [2, 1, 0], 0))
    #print(CG_swig(3, [2, 1, 0], 0, [2, 1, 0], 2, [2, 1, 0], 0))
    
    print(test_CG_swig([1, 0], [2, 0], [3, 0]))
    print(test_CG_swig([2, 0], [3, 0], [3, 0]))
    print(test_CG_swig([2, 0], [3, 0], [5, 0]))

    '''
    A = [1, 0, 0]
    B = [1, 1, 0]
    C = [0, 0, 0]
    D = [1, 1, 0]
    E = [1, 0, 0]
    print(compute_9R(A, B, C, D, E))
    print(math.sqrt(3))
    #debug_CG()
    '''