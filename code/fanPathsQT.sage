"""
    This code manipulates k-tuples (P1, ... Pk) of (m,n)-paths with certain additional conditions:
        - P1, ..., Pk form a fan of non-nesting (m,n)-Dyck paths
        - for all i < j < k, (Pj, Pk) is an interval of the Pi-Tamari lattice
        - we are missing more general conditions that say that for any i1 < ... < ip < j < k, (Pj, Pk) is an interval of a (Pi1, ..., Pip)-Tamari lattice
        
    The goal is that these paths should provide an expression for the character of the multivariate diagonal coinvariant spaces. This is a generalization of the shuffle conjecture.
    Concretely, at the moment, we are given a character H(n,k) by an oracle (Francois Bergeron and Nicolas Thiery), and we want to find a combinatorial formula for this.
"""

""" ========================= Initial data ========================= """

# We need to load the following to make this code work with symmetric and quasisymmetric functions
R.<r,q,t> = QQ['r','q','t']
Sym = SymmetricFunctions(R)
Sym.inject_shorthands()
QSym = QuasiSymmetricFunctions(R)
QSym.inject_shorthands()

Un = s[0]
attach('codeFrancois/Outils_tenseurs.py')
attach('codeFrancois/Rectangular.py')
attach('LLT.py')

from collections import defaultdict 

"""
    Here are the characters given by the oracle.
    They are given by tensors of symmetric functions.
"""

"""
    Note that for (2, n), we could use the following function
    sage: Phi2 = lambda n: add([tensor([s([i]), e([(n+1)//2+i, n//2-i])]) for i in range(n//2+1)])
"""

Phi = {(2, 2): tensor([s([]), s[2]]) + tensor([s[1], s[1,1]]),
       (2, 3): tensor([s([]), s[2, 1]]) + tensor([s[1], s[1, 1, 1]]),
       (2, 4): tensor([s([]), s[2, 2]]) + tensor([s[1], s[2, 1, 1]]) + tensor([s[2], s[1, 1, 1, 1]]),
       (2, 5): tensor([s([]), s[2, 2, 1]]) + tensor([s[1], s[2, 1, 1, 1]]) + tensor([s[2], s[1, 1, 1, 1, 1]]),
       (2, 6): tensor([s([]), s[2, 2, 2]]) + tensor([s[1], s[2, 2, 1, 1]]) + tensor([s[2], s[2, 1, 1, 1, 1]]) + tensor([s[3], s[1, 1, 1, 1, 1, 1]]),
       (3, 2): tensor([s([]), s[2]]) + tensor([s[1], s[1, 1]]),
       (3, 3): tensor([s([]), s[3]]) + tensor([s[1]+s[2], s[2,1]]) + tensor([s[3]+s[1, 1], s[1,1,1]]),
       (3, 4): tensor([s([]), s[3, 1]]) + tensor([s[1], s[2, 1, 1]]) + tensor([s[1], s[2, 2]]) + tensor([s[1, 1], s[1, 1, 1, 1]]) + tensor([s[2], s[2, 1, 1]]) + tensor([s[3], s[1, 1, 1, 1]]),
       (3, 5): tensor([s([]), s[3, 2]]) + tensor([s[1], s[2, 2, 1]]) + tensor([s[1], s[3, 1, 1]]) + tensor([s[1, 1], s[2, 1, 1, 1]]) + tensor([s[2], s[2, 1, 1, 1]]) + tensor([s[2], s[2, 2, 1]]) + tensor([s[2, 1], s[1, 1, 1, 1, 1]]) + tensor([s[3], s[2, 1, 1, 1]]) + tensor([s[4], s[1, 1, 1, 1, 1]]),
       (3, 6): tensor([s([]), s[3, 3]]) + tensor([s[1], s[3, 2, 1]]) + tensor([s[1, 1], s[2, 2, 2]]) + tensor([s[1, 1], s[3, 1, 1, 1]]) + tensor([s[2], s[2, 2, 1, 1]]) + tensor([s[2], s[3, 2, 1]]) + tensor([s[2, 1], s[2, 1, 1, 1, 1]]) + tensor([s[2, 1], s[2, 2, 1, 1]]) + tensor([s[2, 2], s[1, 1, 1, 1, 1, 1]]) + tensor([s[3], s[2, 2, 1, 1]]) + tensor([s[3], s[2, 2, 2]]) + tensor([s[3], s[3, 1, 1, 1]]) + tensor([s[3, 1], s[2, 1, 1, 1, 1]]) + tensor([s[4], s[2, 1, 1, 1, 1]]) + tensor([s[4], s[2, 2, 1, 1]]) + tensor([s[4, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[5], s[2, 1, 1, 1, 1]]) + tensor([s[6], s[1, 1, 1, 1, 1, 1]]),
       (4, 2): tensor([s[1], s[2]]) + tensor([s[2], s[1, 1]]),
       (4, 3): tensor([s([]), s[3]]) + tensor([s[1], s[2, 1]]) + tensor([s[1, 1], s[1, 1, 1]]) + tensor([s[2], s[2, 1]]) + tensor([s[3], s[1, 1, 1]]),
       (4, 4): tensor([s([]), s[4]]) + tensor([s[1]+s[2]+s[3], s[3,1]]) + tensor([s[2]+s[4]+s[2, 1], s[2,2]]) + tensor([s[3]+s[4]+s[5]+s[1, 1]+s[2, 1]+s[3, 1], s[2,1,1]]) + tensor([s[3, 1]+s[4, 1]+s[6]+s[1,1,1], s[1,1,1,1]]),
       (4, 5): tensor([s([]), s[4, 1]]) + tensor([s[1], s[3, 1, 1]]) + tensor([s[1], s[3, 2]]) + tensor([s[1, 1], s[2, 1, 1, 1]]) + tensor([s[1, 1], s[2, 2, 1]]) + tensor([s[1, 1, 1], s[1, 1, 1, 1, 1]]) + tensor([s[2], s[2, 2, 1]]) + tensor([s[2], s[3, 1, 1]]) + tensor([s[2], s[3, 2]]) + tensor([s[2, 1], s[2, 1, 1, 1]]) + tensor([s[2, 1], s[2, 2, 1]]) + tensor([s[3], s[2, 1, 1, 1]]) + tensor([s[3], s[2, 2, 1]]) + tensor([s[3], s[3, 1, 1]]) + tensor([s[3, 1], s[1, 1, 1, 1, 1]]) + tensor([s[3, 1], s[2, 1, 1, 1]]) + tensor([s[4], s[2, 1, 1, 1]]) + tensor([s[4], s[2, 2, 1]]) + tensor([s[4, 1], s[1, 1, 1, 1, 1]]) + tensor([s[5], s[2, 1, 1, 1]]) + tensor([s[6], s[1, 1, 1, 1, 1]]),
       (4, 6): tensor([s[1], s[4, 2]]) + tensor([s[1, 1], s[3, 2, 1]]) + tensor([s[1, 1, 1], s[2, 2, 1, 1]]) + tensor([s[2], s[3, 2, 1]]) + tensor([s[2], s[3, 3]]) + tensor([s[2], s[4, 1, 1]]) + tensor([s[2, 1], s[2, 2, 1, 1]]) + tensor([s[2, 1], s[2, 2, 2]]) + tensor([s[2, 1], s[3, 1, 1, 1]]) + tensor([s[2, 1], s[3, 2, 1]]) + tensor([s[2, 1, 1], s[2, 1, 1, 1, 1]]) + tensor([s[2, 2], s[2, 2, 1, 1]]) + tensor([s[3], s[2, 2, 2]]) + tensor([s[3], s[3, 1, 1, 1]]) + tensor([2*s[3], s[3, 2, 1]]) + tensor([s[3, 1], s[2, 1, 1, 1, 1]]) + tensor([2*s[3, 1], s[2, 2, 1, 1]]) + tensor([s[3, 1], s[2, 2, 2]]) + tensor([s[3, 1], s[3, 1, 1, 1]]) + tensor([s[3, 1, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[3, 2], s[2, 1, 1, 1, 1]]) + tensor([2*s[4], s[2, 2, 1, 1]]) + tensor([s[4], s[3, 1, 1, 1]]) + tensor([s[4], s[3, 2, 1]]) + tensor([2*s[4, 1], s[2, 1, 1, 1, 1]]) + tensor([s[4, 1], s[2, 2, 1, 1]]) + tensor([s[4, 2], s[1, 1, 1, 1, 1, 1]]) + tensor([s[5], s[2, 1, 1, 1, 1]]) + tensor([s[5], s[2, 2, 1, 1]]) + tensor([s[5], s[2, 2, 2]]) + tensor([s[5], s[3, 1, 1, 1]]) + tensor([s[5, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[5, 1], s[2, 1, 1, 1, 1]]) + tensor([s[6], s[2, 1, 1, 1, 1]]) + tensor([s[6], s[2, 2, 1, 1]]) + tensor([s[6, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[7], s[2, 1, 1, 1, 1]]) + tensor([s[8], s[1, 1, 1, 1, 1, 1]]),
       (5, 2): tensor([s[1], s[2]]) + tensor([s[2], s[1, 1]]),
       (5, 3): tensor([s[1], s[3]]) + tensor([s[1, 1], s[2, 1]]) + tensor([s[2], s[2, 1]]) + tensor([s[2, 1], s[1, 1, 1]]) + tensor([s[3], s[2, 1]]) + tensor([s[4], s[1, 1, 1]]),
       (5, 4): tensor([s([]), s[4]]) + tensor([s[1], s[3, 1]]) + tensor([s[1, 1], s[2, 1, 1]]) + tensor([s[1, 1, 1], s[1, 1, 1, 1]]) + tensor([s[2], s[2, 2]]) + tensor([s[2], s[3, 1]]) + tensor([s[2, 1], s[2, 1, 1]]) + tensor([s[2, 1], s[2, 2]]) + tensor([s[3], s[2, 1, 1]]) + tensor([s[3], s[3, 1]]) + tensor([s[3, 1], s[1, 1, 1, 1]]) + tensor([s[3, 1], s[2, 1, 1]]) + tensor([s[4], s[2, 1, 1]]) + tensor([s[4], s[2, 2]]) + tensor([s[4, 1], s[1, 1, 1, 1]]) + tensor([s[5], s[2, 1, 1]]) + tensor([s[6], s[1, 1, 1, 1]]),
       (5, 5): tensor([s([]), s[5]]) + tensor([s[1]+s[2]+s[3]+s[4], s[4,1]]) + tensor([s[2]+s[3]+s[2, 1]+s[4]+s[2, 2]+s[3, 1]+s[5]+s[4, 1]+s[6], s[3,2]]) + tensor([s[1, 1]+s[3]+s[2, 1]+s[4]+2*s[3, 1]+2*s[5]+s[3, 2]+s[4, 1]+s[6]+s[5, 1]+s[7], s[3,1,1]])+ tensor([s[3, 1]+s[6]+2*s[4, 1]+s[4, 2]+s[6, 1]+s[2, 1]+s[4]+s[5]+2*s[5, 1]+s[3, 2]+s[2, 2]+s[7]+s[8]+s[3,1,1]+s[2,1,1], s[2,2,1]]) + tensor([s[3, 1]+s[6]+2*s[4, 1]+s[4, 2]+2*s[6, 1]+s[7, 1]+2*s[5, 1]+s[3, 3]+s[3, 2]+s[7]+s[5, 2]+s[8]+s[9]+s[4,1,1]+s[1,1,1]+s[3,1,1]+s[2,1,1], s[2,1,1,1]]) + tensor([s[1,1,1,1]+s[4,1,1]+s[4,2]+s[3,1,1]+s[6,1]+s[4,3]+s[5,1,1]+s[10]+s[8,1]+s[6,2]+s[7,1], s[1,1,1,1,1]]),
       (5, 6): tensor([s([]), s[5, 1]]) + tensor([s[1], s[4, 1, 1]]) + tensor([s[1], s[4, 2]]) + tensor([s[1, 1], s[3, 1, 1, 1]]) + tensor([s[1, 1], s[3, 2, 1]]) + tensor([s[1, 1, 1], s[2, 1, 1, 1, 1]]) + tensor([s[1, 1, 1], s[2, 2, 1, 1]]) + tensor([s[1, 1, 1, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[2], s[3, 2, 1]]) + tensor([s[2], s[3, 3]]) + tensor([s[2], s[4, 1, 1]]) + tensor([s[2], s[4, 2]]) + tensor([s[2, 1], s[2, 2, 1, 1]]) + tensor([s[2, 1], s[2, 2, 2]]) + tensor([s[2, 1], s[3, 1, 1, 1]]) + tensor([2*s[2, 1], s[3, 2, 1]]) + tensor([s[2, 1], s[3, 3]]) + tensor([s[2, 1, 1], s[2, 1, 1, 1, 1]]) + tensor([s[2, 1, 1], s[2, 2, 1, 1]]) + tensor([s[2, 1, 1], s[2, 2, 2]]) + tensor([s[2, 2], s[2, 2, 1, 1]]) + tensor([s[2, 2], s[3, 2, 1]]) + tensor([s[3], s[3, 1, 1, 1]]) + tensor([2*s[3], s[3, 2, 1]]) + tensor([s[3], s[4, 1, 1]]) + tensor([s[3], s[4, 2]]) + tensor([s[3, 1], s[2, 1, 1, 1, 1]]) + tensor([2*s[3, 1], s[2, 2, 1, 1]]) + tensor([s[3, 1], s[2, 2, 2]]) + tensor([2*s[3, 1], s[3, 1, 1, 1]]) + tensor([2*s[3, 1], s[3, 2, 1]]) + tensor([s[3, 1, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[3, 1, 1], s[2, 1, 1, 1, 1]]) + tensor([s[3, 1, 1], s[2, 2, 1, 1]]) + tensor([s[3, 2], s[2, 1, 1, 1, 1]]) + tensor([s[3, 2], s[2, 2, 1, 1]]) + tensor([s[3, 2], s[2, 2, 2]]) + tensor([s[3, 2], s[3, 1, 1, 1]]) + tensor([s[3, 3], s[2, 1, 1, 1, 1]]) + tensor([s[4], s[2, 2, 1, 1]]) + tensor([s[4], s[2, 2, 2]]) + tensor([s[4], s[3, 1, 1, 1]]) + tensor([2*s[4], s[3, 2, 1]]) + tensor([s[4], s[3, 3]]) + tensor([s[4], s[4, 1, 1]]) + tensor([2*s[4, 1], s[2, 1, 1, 1, 1]]) + tensor([3*s[4, 1], s[2, 2, 1, 1]]) + tensor([s[4, 1], s[2, 2, 2]]) + tensor([s[4, 1], s[3, 1, 1, 1]]) + tensor([s[4, 1], s[3, 2, 1]]) + tensor([s[4, 1, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[4, 1, 1], s[2, 1, 1, 1, 1]]) + tensor([s[4, 2], s[1, 1, 1, 1, 1, 1]]) + tensor([s[4, 2], s[2, 1, 1, 1, 1]]) + tensor([s[4, 2], s[2, 2, 1, 1]]) + tensor([s[4, 3], s[1, 1, 1, 1, 1, 1]]) + tensor([s[5], s[2, 2, 1, 1]]) + tensor([s[5], s[2, 2, 2]]) + tensor([2*s[5], s[3, 1, 1, 1]]) + tensor([2*s[5], s[3, 2, 1]]) + tensor([2*s[5, 1], s[2, 1, 1, 1, 1]]) + tensor([2*s[5, 1], s[2, 2, 1, 1]]) + tensor([s[5, 1], s[2, 2, 2]]) + tensor([s[5, 1], s[3, 1, 1, 1]]) + tensor([s[5, 1, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[5, 2], s[2, 1, 1, 1, 1]]) + tensor([s[6], s[2, 1, 1, 1, 1]]) + tensor([2*s[6], s[2, 2, 1, 1]]) + tensor([s[6], s[3, 1, 1, 1]]) + tensor([s[6], s[3, 2, 1]]) + tensor([s[6, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([2*s[6, 1], s[2, 1, 1, 1, 1]]) + tensor([s[6, 1], s[2, 2, 1, 1]]) + tensor([s[6, 2], s[1, 1, 1, 1, 1, 1]]) + tensor([s[7], s[2, 1, 1, 1, 1]]) + tensor([s[7], s[2, 2, 1, 1]]) + tensor([s[7], s[2, 2, 2]]) + tensor([s[7], s[3, 1, 1, 1]]) + tensor([s[7, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[7, 1], s[2, 1, 1, 1, 1]]) + tensor([s[8], s[2, 1, 1, 1, 1]]) + tensor([s[8], s[2, 2, 1, 1]]) + tensor([s[8, 1], s[1, 1, 1, 1, 1, 1]]) + tensor([s[9], s[2, 1, 1, 1, 1]]) + tensor([s[10], s[1, 1, 1, 1, 1, 1]]),
       (6, 2): tensor([s[2], s[2]]) + tensor([s[3], s[1, 1]]),
       (6, 3): tensor([s[1, 1], s[3]]) + tensor([s[2, 1], s[2, 1]]) + tensor([s[2, 2], s[1, 1, 1]]) + tensor([s[3], s[3]]) + tensor([s[3, 1], s[2, 1]]) + tensor([s[4], s[2, 1]]) + tensor([s[4, 1], s[1, 1, 1]]) + tensor([s[5], s[2, 1]]) + tensor([s[6], s[1, 1, 1]]),
       (6, 4): tensor([s[1, 1, 1], s[2, 2]]) + tensor([s[2], s[4]]) + tensor([s[2, 1], s[3, 1]]) + tensor([s[2, 1, 1], s[2, 1, 1]]) + tensor([s[2, 2], s[2, 2]]) + tensor([s[3], s[3, 1]]) + tensor([s[3, 1], s[2, 1, 1]]) + tensor([s[3, 1], s[2, 2]]) + tensor([s[3, 1], s[3, 1]]) + tensor([s[3, 1, 1], s[1, 1, 1, 1]]) + tensor([s[3, 2], s[2, 1, 1]]) + tensor([s[4], s[2, 2]]) + tensor([s[4], s[3, 1]]) + tensor([2*s[4, 1], s[2, 1, 1]]) + tensor([s[4, 1], s[2, 2]]) + tensor([s[4, 2], s[1, 1, 1, 1]]) + tensor([s[5], s[2, 1, 1]]) + tensor([s[5], s[3, 1]]) + tensor([s[5, 1], s[1, 1, 1, 1]]) + tensor([s[5, 1], s[2, 1, 1]]) + tensor([s[6], s[2, 1, 1]]) + tensor([s[6], s[2, 2]]) + tensor([s[6, 1], s[1, 1, 1, 1]]) + tensor([s[7], s[2, 1, 1]]) + tensor([s[8], s[1, 1, 1, 1]]),
       (6, 5): tensor([s([]), s[5]]) + tensor([s[1], s[4, 1]]) + tensor([s[1, 1], s[3, 1, 1]]) + tensor([s[1, 1, 1], s[2, 1, 1, 1]]) + tensor([s[1, 1, 1, 1], s[1, 1, 1, 1, 1]]) + tensor([s[2], s[3, 2]]) + tensor([s[2], s[4, 1]]) + tensor([s[2, 1], s[2, 2, 1]]) + tensor([s[2, 1], s[3, 1, 1]]) + tensor([s[2, 1], s[3, 2]]) + tensor([s[2, 1, 1], s[2, 1, 1, 1]]) + tensor([s[2, 1, 1], s[2, 2, 1]]) + tensor([s[2, 2], s[2, 2, 1]]) + tensor([s[2, 2], s[3, 2]]) + tensor([s[3], s[3, 1, 1]]) + tensor([s[3], s[3, 2]]) + tensor([s[3], s[4, 1]]) + tensor([s[3, 1], s[2, 1, 1, 1]]) + tensor([s[3, 1], s[2, 2, 1]]) + tensor([2*s[3, 1], s[3, 1, 1]]) + tensor([s[3, 1], s[3, 2]]) + tensor([s[3, 1, 1], s[1, 1, 1, 1, 1]]) + tensor([s[3, 1, 1], s[2, 1, 1, 1]]) + tensor([s[3, 1, 1], s[2, 2, 1]]) + tensor([s[3, 2], s[2, 1, 1, 1]]) + tensor([s[3, 2], s[2, 2, 1]]) + tensor([s[3, 2], s[3, 1, 1]]) + tensor([s[3, 3], s[2, 1, 1, 1]]) + tensor([s[4], s[2, 2, 1]]) + tensor([s[4], s[3, 1, 1]]) + tensor([s[4], s[3, 2]]) + tensor([s[4], s[4, 1]]) + tensor([2*s[4, 1], s[2, 1, 1, 1]]) + tensor([2*s[4, 1], s[2, 2, 1]]) + tensor([s[4, 1], s[3, 1, 1]]) + tensor([s[4, 1], s[3, 2]]) + tensor([s[4, 1, 1], s[1, 1, 1, 1, 1]]) + tensor([s[4, 1, 1], s[2, 1, 1, 1]]) + tensor([s[4, 2], s[1, 1, 1, 1, 1]]) + tensor([s[4, 2], s[2, 1, 1, 1]]) + tensor([s[4, 2], s[2, 2, 1]]) + tensor([s[4, 3], s[1, 1, 1, 1, 1]]) + tensor([s[5], s[2, 2, 1]]) + tensor([2*s[5], s[3, 1, 1]]) + tensor([s[5], s[3, 2]]) + tensor([2*s[5, 1], s[2, 1, 1, 1]]) + tensor([2*s[5, 1], s[2, 2, 1]]) + tensor([s[5, 1], s[3, 1, 1]]) + tensor([s[5, 1, 1], s[1, 1, 1, 1, 1]]) + tensor([s[5, 2], s[2, 1, 1, 1]]) + tensor([s[6], s[2, 1, 1, 1]]) + tensor([s[6], s[2, 2, 1]]) + tensor([s[6], s[3, 1, 1]]) + tensor([s[6], s[3, 2]]) + tensor([s[6, 1], s[1, 1, 1, 1, 1]]) + tensor([2*s[6, 1], s[2, 1, 1, 1]]) + tensor([s[6, 1], s[2, 2, 1]]) + tensor([s[6, 2], s[1, 1, 1, 1, 1]]) + tensor([s[7], s[2, 1, 1, 1]]) + tensor([s[7], s[2, 2, 1]]) + tensor([s[7], s[3, 1, 1]]) + tensor([s[7, 1], s[1, 1, 1, 1, 1]]) + tensor([s[7, 1], s[2, 1, 1, 1]]) + tensor([s[8], s[2, 1, 1, 1]]) + tensor([s[8], s[2, 2, 1]]) + tensor([s[8, 1], s[1, 1, 1, 1, 1]]) + tensor([s[9], s[2, 1, 1, 1]]) + tensor([s[10], s[1, 1, 1, 1, 1]]),
       (6, 6): tensor([s([]), s[6]]) + tensor([s[1]+s[2]+s[3]+s[4]+s[5], s[5,1]]) + tensor([s[2]+s[2,1]+s[2,2]+s[3]+s[3,1]+s[3,2]+2*s[4]+2*s[4,1]+s[4,2]+s[5]+s[5,1]+2*s[6]+s[6,1]+s[7]+s[8], s[4,2]]) + tensor([s[1,1]+s[2,1]+s[3]+2*s[3,1]+s[3,2]+s[3,3]+s[4]+2*s[4,1]+s[4,2]+2*s[5]+2*s[5,1]+s[5,2]+2*s[6]+s[6,1]+2*s[7]+s[7,1]+s[8]+s[9], s[4,1,1]]) + tensor([s[2,2]+s[2,2,1]+s[3]+s[3,1]+s[3,2]+s[4,1]+s[4,1,1]+s[4,2]+s[5]+s[5,1]+s[5,2]+s[6]+s[6,1]+s[7]+s[7,1]+s[9], s[3,3]]) + tensor([s[2,2]+s[2,2,1]+s[3,1,1]+s[3,1,1,1]+s[3,2,1]+s[3,3,1]+s[4,1]+s[4,1,1]+2*s[4,2]+s[4,2,1]+s[4,3]+s[4,4]+s[5,1]+2*s[5,1,1]+2*s[5,2]+s[5,2,1]+s[5,3]+s[6]+2*s[6,1]+s[6,1,1]+2*s[6,2]+s[6,3]+2*s[7,1]+s[7,1,1]+s[7,2]+s[8]+2*s[8,1]+s[8,2]+s[9]+s[9,1]+s[10]+s[10,1]+s[12], s[2,2,2]]) + tensor([s[1,1,1]+s[2,1,1]+s[3,1]+2*s[3,1,1]+s[3,2]+s[3,2,1]+2*s[3,3]+s[3,3,1]+2*s[4,1]+2*s[4,1,1]+2*s[4,2]+s[4,2,1]+2*s[4,3]+3*s[5,1]+2*s[5,1,1]+3*s[5,2]+s[5,2,1]+2*s[5,3]+s[6]+4*s[6,1]+s[6,1,1]+3*s[6,2]+s[6,3]+s[7]+4*s[7,1]+s[7,1,1]+2*s[7,2]+2*s[8]+3*s[8,1]+s[8,2]+2*s[9]+2*s[9,1]+2*s[10]+s[10,1]+s[11]+s[12], s[3,1,1,1]]) + tensor([s[2,1]+s[2,1,1]+s[2,2]+s[2,2,1]+2*s[3,1]+2*s[3,1,1]+3*s[3,2]+2*s[3,2,1]+s[3,3]+s[4]+3*s[4,1]+2*s[4,1,1]+4*s[4,2]+s[4,2,1]+2*s[4,3]+2*s[5]+5*s[5,1]+2*s[5,1,1]+4*s[5,2]+s[5,3]+2*s[6]+5*s[6,1]+s[6,1,1]+3*s[6,2]+3*s[7]+4*s[7,1]+s[7,2]+3*s[8]+3*s[8,1]+2*s[9]+s[9,1]+2*s[10]+s[11], s[3,2,1]]) + tensor([s[1,1,1,1]+s[2,1,1,1]+s[3,1,1]+s[3,1,1,1]+s[3,2,1]+s[3,3,1]+2*s[4,1,1]+s[4,1,1,1]+s[4,2]+2*s[4,2,1]+2*s[4,3]+2*s[4,3,1]+s[4,4]+3*s[5,1,1]+s[5,1,1,1]+s[5,2]+2*s[5,2,1]+2*s[5,3]+s[5,3,1]+s[5,4]+s[6,1]+3*s[6,1,1]+3*s[6,2]+2*s[6,2,1]+3*s[6,3]+s[6,4]+2*s[7,1]+3*s[7,1,1]+3*s[7,2]+s[7,2,1]+2*s[7,3]+3*s[8,1]+2*s[8,1,1]+3*s[8,2]+s[8,3]+3*s[9,1]+s[9,1,1]+2*s[9,2]+s[10]+3*s[10,1]+s[10,2]+s[11]+2*s[11,1]+s[12]+s[12,1]+s[13]+s[14], s[2,1,1,1,1]]) + tensor([s[2,1,1]+s[2,1,1,1]+s[2,2,1]+s[3,1,1]+s[3,1,1,1]+s[3,2]+2*s[3,2,1]+s[3,3]+s[3,3,1]+s[4,1]+3*s[4,1,1]+s[4,1,1,1]+2*s[4,2]+3*s[4,2,1]+2*s[4,3]+s[4,3,1]+s[4,4]+2*s[5,1]+3*s[5,1,1]+4*s[5,2]+2*s[5,2,1]+3*s[5,3]+s[5,4]+3*s[6,1]+4*s[6,1,1]+4*s[6,2]+s[6,2,1]+2*s[6,3]+s[7]+4*s[7,1]+2*s[7,1,1]+4*s[7,2]+s[7,3]+s[8]+4*s[8,1]+s[8,1,1]+2*s[8,2]+2*s[9]+4*s[9,1]+s[9,2]+s[10]+2*s[10,1]+2*s[11]+s[11,1]+s[12]+s[13], s[2,2,1,1]]) + tensor([2*s[8,1,1]+s[8,2,1]+s[8,2]+s[6,2,1]+s[9,1,1]+s[15]+s[9,3]+s[6,1,1]+s[13,1]+s[10,2]+s[4,1,1,1]+s[6,1,1,1]+s[7,2,1]+s[8,3]+s[12,1]+s[4,3,1]+s[4,4]+s[4,2,1]+s[9,2]+s[11,2]+s[11,1]+s[4,4,1]+s[6,4]+s[5,1,1,1]+s[7,2]+s[5,3,1]+s[7,3]+s[10,1]+s[7,1,1]+s[6,3]+s[1,1,1,1,1]+s[3,1,1,1]+s[7,4]+s[10,1,1]+s[5,2,1]+s[6,3,1], s[1,1,1,1,1,1]]),
       (7, 4): tensor([s[1, 1], s[4]]) + tensor([s[1, 1, 1], s[3, 1]]) + tensor([s[2, 1], s[3, 1]]) + tensor([s[2, 1, 1], s[2, 1, 1]]) + tensor([s[2, 1, 1], s[2, 2]]) + tensor([s[2, 2], s[2, 1, 1]]) + tensor([s[2, 2, 1], s[1, 1, 1, 1]]) + tensor([s[3], s[4]]) + tensor([s[3, 1], s[2, 2]]) + tensor([2*s[3, 1], s[3, 1]]) + tensor([s[3, 1, 1], s[2, 1, 1]]) + tensor([s[3, 2], s[2, 1, 1]]) + tensor([s[3, 2], s[2, 2]]) + tensor([s[4], s[3, 1]]) + tensor([2*s[4, 1], s[2, 1, 1]]) + tensor([s[4, 1], s[2, 2]]) + tensor([s[4, 1], s[3, 1]]) + tensor([s[4, 1, 1], s[1, 1, 1, 1]]) + tensor([s[4, 2], s[1, 1, 1, 1]]) + tensor([s[4, 2], s[2, 1, 1]]) + tensor([s[5], s[2, 2]]) + tensor([s[5], s[3, 1]]) + tensor([2*s[5, 1], s[2, 1, 1]]) + tensor([s[5, 1], s[2, 2]]) + tensor([s[5, 2], s[1, 1, 1, 1]]) + tensor([s[6], s[2, 1, 1]]) + tensor([s[6], s[3, 1]]) + tensor([s[6, 1], s[1, 1, 1, 1]]) + tensor([s[6, 1], s[2, 1, 1]]) + tensor([s[7], s[2, 1, 1]]) + tensor([s[7], s[2, 2]]) + tensor([s[7, 1], s[1, 1, 1, 1]]) + tensor([s[8], s[2, 1, 1]]) + tensor([s[9], s[1, 1, 1, 1]]),
       (8, 4): tensor([s[1, 1, 1], s[4]]) + tensor([s[2, 1, 1], s[3, 1]]) + tensor([s[2, 2, 1], s[2, 1, 1]]) + tensor([s[2, 2, 2], s[1, 1, 1, 1]]) + tensor([s[3, 1], s[4]]) + tensor([s[3, 1, 1], s[2, 2]]) + tensor([s[3, 1, 1], s[3, 1]]) + tensor([s[3, 2], s[3, 1]]) + tensor([s[3, 2, 1], s[2, 1, 1]]) + tensor([s[3, 2, 1], s[2, 2]]) + tensor([s[3, 3], s[3, 1]]) + tensor([s[4, 1], s[3, 1]]) + tensor([s[4, 1], s[4]]) + tensor([s[4, 1, 1], s[2, 1, 1]]) + tensor([s[4, 1, 1], s[3, 1]]) + tensor([s[4, 2], s[2, 1, 1]]) + tensor([s[4, 2], s[2, 2]]) + tensor([s[4, 2], s[3, 1]]) + tensor([s[4, 2, 1], s[1, 1, 1, 1]]) + tensor([s[4, 2, 1], s[2, 1, 1]]) + tensor([s[4, 3], s[2, 1, 1]]) + tensor([s[4, 3], s[2, 2]]) + tensor([s[4, 4], s[1, 1, 1, 1]]) + tensor([s[5, 1], s[2, 2]]) + tensor([2*s[5, 1], s[3, 1]]) + tensor([s[5, 1, 1], s[2, 1, 1]]) + tensor([s[5, 1, 1], s[2, 2]]) + tensor([2*s[5, 2], s[2, 1, 1]]) + tensor([s[5, 2], s[2, 2]]) + tensor([s[5, 2], s[3, 1]]) + tensor([s[5, 2, 1], s[1, 1, 1, 1]]) + tensor([s[5, 3], s[2, 1, 1]]) + tensor([s[6], s[4]]) + tensor([s[6, 1], s[2, 1, 1]]) + tensor([s[6, 1], s[2, 2]]) + tensor([2*s[6, 1], s[3, 1]]) + tensor([s[6, 1, 1], s[2, 1, 1]]) + tensor([s[6, 2], s[1, 1, 1, 1]]) + tensor([2*s[6, 2], s[2, 1, 1]]) + tensor([s[6, 2], s[2, 2]]) + tensor([s[6, 3], s[1, 1, 1, 1]]) + tensor([s[7], s[3, 1]]) + tensor([2*s[7, 1], s[2, 1, 1]]) + tensor([s[7, 1], s[2, 2]]) + tensor([s[7, 1], s[3, 1]]) + tensor([s[7, 1, 1], s[1, 1, 1, 1]]) + tensor([s[7, 2], s[1, 1, 1, 1]]) + tensor([s[7, 2], s[2, 1, 1]]) + tensor([s[8], s[2, 2]]) + tensor([s[8], s[3, 1]]) + tensor([2*s[8, 1], s[2, 1, 1]]) + tensor([s[8, 1], s[2, 2]]) + tensor([s[8, 2], s[1, 1, 1, 1]]) + tensor([s[9], s[2, 1, 1]]) + tensor([s[9], s[3, 1]]) + tensor([s[9, 1], s[1, 1, 1, 1]]) + tensor([s[9, 1], s[2, 1, 1]]) + tensor([s[10], s[2, 1, 1]]) + tensor([s[10], s[2, 2]]) + tensor([s[10, 1], s[1, 1, 1, 1]]) + tensor([s[11], s[2, 1, 1]]) + tensor([s[12], s[1, 1, 1, 1]]),
       (9, 4): tensor([s[1, 1, 1], s[4]]) + tensor([s[2, 1, 1], s[3, 1]]) + tensor([s[2, 2, 1], s[2, 1, 1]]) + tensor([s[2, 2, 2], s[1, 1, 1, 1]]) + tensor([s[3, 1], s[4]]) + tensor([s[3, 1, 1], s[2, 2]]) + tensor([s[3, 1, 1], s[3, 1]]) + tensor([s[3, 2], s[3, 1]]) + tensor([s[3, 2, 1], s[2, 1, 1]]) + tensor([s[3, 2, 1], s[2, 2]]) + tensor([s[3, 3], s[3, 1]]) + tensor([s[4, 1], s[3, 1]]) + tensor([s[4, 1], s[4]]) + tensor([s[4, 1, 1], s[2, 1, 1]]) + tensor([s[4, 1, 1], s[3, 1]]) + tensor([s[4, 2], s[2, 1, 1]]) + tensor([s[4, 2], s[2, 2]]) + tensor([s[4, 2], s[3, 1]]) + tensor([s[4, 2, 1], s[1, 1, 1, 1]]) + tensor([s[4, 2, 1], s[2, 1, 1]]) + tensor([s[4, 3], s[2, 1, 1]]) + tensor([s[4, 3], s[2, 2]]) + tensor([s[4, 4], s[1, 1, 1, 1]]) + tensor([s[5, 1], s[2, 2]]) + tensor([2*s[5, 1], s[3, 1]]) + tensor([s[5, 1, 1], s[2, 1, 1]]) + tensor([s[5, 1, 1], s[2, 2]]) + tensor([2*s[5, 2], s[2, 1, 1]]) + tensor([s[5, 2], s[2, 2]]) + tensor([s[5, 2], s[3, 1]]) + tensor([s[5, 2, 1], s[1, 1, 1, 1]]) + tensor([s[5, 3], s[2, 1, 1]]) + tensor([s[6], s[4]]) + tensor([s[6, 1], s[2, 1, 1]]) + tensor([s[6, 1], s[2, 2]]) + tensor([2*s[6, 1], s[3, 1]]) + tensor([s[6, 1, 1], s[2, 1, 1]]) + tensor([s[6, 2], s[1, 1, 1, 1]]) + tensor([2*s[6, 2], s[2, 1, 1]]) + tensor([s[6, 2], s[2, 2]]) + tensor([s[6, 3], s[1, 1, 1, 1]]) + tensor([s[7], s[3, 1]]) + tensor([2*s[7, 1], s[2, 1, 1]]) + tensor([s[7, 1], s[2, 2]]) + tensor([s[7, 1], s[3, 1]]) + tensor([s[7, 1, 1], s[1, 1, 1, 1]]) + tensor([s[7, 2], s[1, 1, 1, 1]]) + tensor([s[7, 2], s[2, 1, 1]]) + tensor([s[8], s[2, 2]]) + tensor([s[8], s[3, 1]]) + tensor([2*s[8, 1], s[2, 1, 1]]) + tensor([s[8, 1], s[2, 2]]) + tensor([s[8, 2], s[1, 1, 1, 1]]) + tensor([s[9], s[2, 1, 1]]) + tensor([s[9], s[3, 1]]) + tensor([s[9, 1], s[1, 1, 1, 1]]) + tensor([s[9, 1], s[2, 1, 1]]) + tensor([s[10], s[2, 1, 1]]) + tensor([s[10], s[2, 2]]) + tensor([s[10, 1], s[1, 1, 1, 1]]) + tensor([s[11], s[2, 1, 1]]) + tensor([s[12], s[1, 1, 1, 1]])}

"""
To evaluate these tensors, we just need to use plethism...
"""

""" ========================= Binomials ========================= """

# Return the abstract binomial binom(var, j).
def abstractBinomial(var, j):
    # SR is the symbolic ring. The option hold=True avoids the evaluation of the binomial coefficient.
    return SR(var).binomial(j, hold=True)

# Change the polynomial basis r^j to the binomial basis binom(r,j).
# More precisely, return the coefficients of P(r) in the binomial basis binomial(r+k, i).
def expressPolyOnBinomialCoeffs(P, k=0):
    coefficients = []
    for i in range(P.degree(),-1,-1):
        coeff = P.coefficient(r^i) * factorial(i)
        coefficients.append(coeff)
        P -= coeff * binomial(r+k, i)
    coefficients.reverse()
    return coefficients

def expressSymFunctOnBinomialCoeffs(f, k=0):
    return [(mu, expressPolyOnBinomialCoeffs(coeff, k=k)) for (mu, coeff) in f]

""" ========================= Chains in a poset ========================= """

# Recall that the strict chains of a poset P are given by P.chains().
# If you want the strict chains of length i, you use P.chains().elements_of_depth_iterator(i).

# Return the weak chains of length k in the poset P.
# For that, we enumerate all strict chains of length at most k and repeat some of their elements.
def weakChains(P, k, element_constructor=list):
    weakChains = []
    for i in range(1,k+1):
        for strictChain in P.chains().elements_of_depth_iterator(i):
            for composition in Compositions(k, length=i):
                weakChains.append(element_constructor(flatten([[strictChain[j]]*composition[j] for j in range(i)])))
    return weakChains

# Compute the strict chains in a poset P, with endpoints bp and tp.
def strictChainsInterval(P, bp, tp, element_constructor=tuple):
    return P.subposet(P.interval(bp,tp)[1:-1]).chains(element_constructor=element_constructor)

# Compute the strict chains in a poset P in the intervals between chain[i] and chain[i+1] for all i.
def strictChainsIntervals(P, chain, element_constructor=tuple):
    if len(chain) == 0:
        return []
    # compute all chains between chain[i] and chain[i+1] for all i
    return [strictChainsInterval(P, chain[i], chain[i+1], element_constructor=element_constructor) for i in range(len(chain)-1)]

# Construct a superChain from an initial chain and some chains in the intervals between chain[i] and chain[i+1] for all i.
def constructChain(chain, chainsIntervals):
    return flatten(zip(map(lambda x: tuple([x]), chain), chainsIntervals)) + [chain[-1]]

# Compute the strict superchains of a given chain in a poset P, with same endpoints as the initial chain.
# This is just the Cartesian product of the chains between two consecutive elements of the chain.
def strictBoundedSuperChains(P, chain, element_constructor=list):
    if len(chain) == 0:
        return P.chains()
    # compute all chains between chain[i] and chain[i+1] for all i
    # not included if chain[i] == chain[i+1]
    chainsIntervals = [[(chain[i],) + c for c in strictChainsInterval(P, chain[i], chain[i+1], element_constructor=tuple)] for i in range(len(chain)-1) if chain[i] != chain[i+1]]
    return [element_constructor(flatten(x)) for x in cartesian_product(chainsIntervals + [[chain[-1]]])]

# Compute the strict superchains of a given chain in a poset P.
# For that, we compute:
#    - the lower ideal of the bottom element of chain,
#    - the intervals between two consecutive elements in chain,
#    - the upper ideal of the top element of chain.
# The set of superchains of chain is the cartesian product of all chains in the resulting posets.
def strictSuperChains(P, chain, element_constructor=list):
    if len(chain) == 0:
        return P.chains()
    return [element_constructor(flatten(x)) for x in cartesian_product([P.subposet(P.principal_lower_set(chain[0])[:-1]).chains()] + [strictBoundedSuperChains(P, chain, element_constructor=tuple)] + [P.subposet(P.principal_upper_set(chain[-1])[1:]).chains()])]

""" ========================= Maxima and maximal elements in a list ========================= """

# Test whether a list l1 is a is_sublist of a list l2.
def is_sublist(l1, l2):
    return (len(l1) == 0) or (len(l2) != 0 and ((l1[0] == l2[0] and is_sublist(l1[1:], l2[1:])) or is_sublist(l1, l2[1:])))

# Record the maximum of a list and the number of occurences of such maxima.
def recordMaxima(l):
    max = l[0]
    nMax = 0
    for x in l:
        if x == max:
            nMax += 1
        if x > max:
            max = x
            nMax = 1
    return (max, nMax)

# Keep the maximal elements of a list with respect to a functional f.
def keepMaximalElements(l,f):
    maxValue = f(l[0])
    maxElements = [l[0]]
    for x in l:
        if f(x) == maxValue:
            maxElements.append(x)
        if f(x) > maxValue:
            maxValue = f(x)
            maxElements = [x]
    return (maxValue, maxElements)

""" ========================= Manipulations of paths and bracket vectors ========================= """

# Reverse the path (mirror and change 1 to 0).
def reversePath(path):
    return [1-path[len(path)-1-i] for i in range(len(path))]

# Return the partition of the path (which is visible above the path).
def partitionPath(path):
    partition = [len(path)//2]
    for i in range(len(path)):
        if path[len(path) - 1 -i] == 0:
            partition[-1] -= 1
        else:
            partition.append(partition[-1])
    return tuple(partition)

# Return the path of a partition (the limit of the partition).
def pathPartition(partition):
    n = len(partition)
    enlargedPartition = list(reversed(partition)) + [n]
    return Word(flatten([[1]+[0]*(enlargedPartition[i+1] - enlargedPartition[i]) for i in range(n)]))

# Return the sequence of the lengths of the up steps of the path.
# That is, if path = N^h1 E^w1 N^h2 E^w2 N^h3 ..., the up steps sequence is (h1, h2, ...).
def upStepsVector(path):
    upSteps = []
    height = 0
    for step in path:
        if step == 1:
            height += 1
        else:
            if height != 0:
                upSteps.append(height)
            height = 0
    return upSteps

@cached_function
# Return the partition of the lengths of the up steps of the path.
def upStepsPartition(path):
    return Partition(reversed(sorted(upStepsVector(path))))

@cached_function
# Return the dictionnary of (m,n)-Dyck paths according to their upStepsPartition.
def upStepDictionnary(m,n):
    upStepDictionnary = {}
    for path in TamariLattice(m,n):
        partition = upStepsPartition(path)
        if not upStepDictionnary.has_key(partition):
            upStepDictionnary[partition] = []
        upStepDictionnary[partition].append(path)
    return upStepDictionnary

# Return the sequence of the lengths of the east steps of the path.
# That is, if path = N^h1 E^w1 N^h2 E^w2 N^h3 ..., the east steps sequence is (wn, ..., w2, w1).
def eastStepsVector(path):
    eastSteps = []
    width = 0
    for i in range(len(path)):
        if path[len(path)-1-i] == 0:
            width += 1
        else:
            if width != 0:
                eastSteps.append(width)
            width = 0
    return eastSteps

@cached_function
# Return the partition of the lengths of the east steps of the path.
def eastStepsPartition(path):
    return Partition(reversed(sorted(eastStepsVector(path))))

# Return the width sequence of the path.
# That is, if path = N E^w1 N E^w2 N ..., the width sequence is (w1, w2, ...).
def widths(path):
    widths = [0]
    for step in path:
        if step == 1:
            widths.append(0)
        else:
            widths[-1] += 1
    return widths

# Return the horizontal bracket vector of a path path2 with respect to the path path1.
def bracketVectorWidth(path1, path2):
    widths1 = widths(path1)
    widths2 = widths(path2)
    bracketVectorWidth = [-1]*(len(path1)+1)
    topPosition = -1
    for i in range(len(widths1)):
        topPosition += widths1[i]
        currentPosition = topPosition
        j = 0
        while j < widths2[i]:
            if bracketVectorWidth[currentPosition] == -1:
                bracketVectorWidth[currentPosition] = i
                j += 1
            currentPosition -= 1
    return bracketVectorWidth

# Return the height sequence of the path.
# That is, if path = N^h1 E N^h2 E ..., the heights sequence is (h1, h2, ...).
def heights(path):
    heights = [0]
    for step in reversed(path):
        if step == 0:
            heights.append(0)
        else:
            heights[-1] += 1
    return heights

# Return the vertical bracket vector of a path path2 with respect to the path path1.
def bracketVectorHeight(path1, path2):
    heights1 = heights(path1)
    heights2 = heights(path2)
    bracketVectorHeight = [-1]*(len(path1)+1)
    topPosition = -1
    for i in range(len(heights1)):
        topPosition += heights1[i]
        currentPosition = topPosition
        j = 0
        while j < heights2[i]:
            if bracketVectorHeight[currentPosition] == -1:
                bracketVectorHeight[currentPosition] = i
                j += 1
            currentPosition -= 1
    return bracketVectorHeight

# Return the sequence of (x,y) positions of the enpoints of the steps of the path
def positions(path):
    currentPosition = [0,0]
    positions = [[0,0]]
    for step in path:
        currentPosition[step] += 1
        positions.append(tuple(currentPosition))
    return positions

# Compare two vectors v1 and v2 componentwise.
def componentwiseLess(v1, v2):
    return all([v1[i] <= v2[i] for i in range(len(v1))])

# Return the words with i zeros and j ones.
def wordSubsets(i,j):
    return [Word([int(k in s) for k in range(i+j)]) for s in Subsets(range(i+j),j)]

# zeta map
# sends the (area,dinv) to (dinv,bounce)
def zeta(m, n, path):
    shm = int(gcd(m,n) != 1)
    path = path + [0]*shm
    return map(lambda x: x[2], sorted(zip(map(lambda (x,y):y*(m+shm)-x*n, positions(path)[:-1]), range(len(path)), path)))[:len(path)+(-1)*shm]

# eta map
# sends the (area,dinv) to (dinv,bounce)
def eta(m, n, path):
    shm = int(gcd(m,n) != 1)
    path = path + [0]*shm
    return list(reversed(map(lambda x: x[2], sorted(zip(map(lambda (x,y):y*(m+shm)-x*n, positions(path)[1:]), range(len(path)), path)))))[:len(path)+(-1)*shm]

""" ========================= drawings ========================= """

# Represent a path in the terminal
def latticeRepresentation(path):
    repr = ""
    dist = 0
    for i in range(len(path)):
        if path[i] == 1:
            repr +=  "\n" + " "*dist + "|"
        else:
            repr += "_"
            dist += 1
    return repr

# Draw a chain of paths in the plane.
def draw(words, list_colors=[], thickness=2, gridsize=(1,1), shifts=(.1,.1)):
    # import plot facilities
    from sage.plot.line import line
    from sage.plot.colors import colors
    # prepare colors
    list_colors += ['red', 'blue', 'green', 'orange', 'yellow', 'purple']
    list_colors += list(colors)
    # prepare lines
    paths = line([(1, 1)])
    for i in range(len(words)):
        (x,y) = (-i*shifts[0]*gridsize[0], i*shifts[1]*gridsize[1])
        path = [(x,y)]
        for j in range(len(words[i])):
            x += (1-words[i][j])*gridsize[0]
            y += words[i][j] * gridsize[1]
            path.append((x,y))
        paths += line(path, color=list_colors[i], thickness=thickness)
    paths.axes(False)
    paths.set_aspect_ratio(1)
    return paths

# Draw a chain of paths after reversing.
def drawWordsTranspose(words, **args):
    return draw([reversePath([int(c) for c in list(word)]) for word in words], **args)

""" ========================= nu-Tamari lattice ========================= """

@cached_function
# Return the top path in the (m,n)-Tamari lattice.
def topPath(m,n):
    return Word([1]*n + [0]*m)

@cached_function
# Return the bottom path in the (m,n)-Tamari lattice.
def bottomPath(m,n):
    res = []
    j = 0
    for i in range(1,n+1):
        res += [1] + [0]*(floor(i*m/n)-j)
        j = floor(i*m/n)
    return Word(res)

# Is p smaller than q in nu-Tamari?
def compareNuTamari(nu, p, q):
    return componentwiseLess(bracketVectorHeight(nu, p), bracketVectorHeight(nu, q))

# Return the nu-Tamari lattice.
def nuTamariLattice(nu):
    elem = [q for q in wordSubsets(len(nu)-sum(nu), sum(nu)) if compareNuTamari(nu, nu, q)]
    rel = [[p,q] for p in elem for q in elem if compareNuTamari(nu, p, q)]
    return Poset([elem, rel])

@cached_function
# Return the (m,n)-Tamari lattice.
def TamariLattice(m,n):
    return nuTamariLattice(bottomPath(m,n))

""" ========================= Stanley poset ========================= """

# Is p below q? (path inclusion)
def isBelow(p,q):
    return componentwiseLess(partitionPath(q), partitionPath(p))

# Return the Stanley poset, that is, the poset of Dyck path ordered by inclusion.
def StanleyPoset(m, n):
    elem = wordSubsets(m, n)
    rel = [[p,q] for p in elem for q in elem if isBelow(p,q)]
    return Poset([elem, rel])

""" ========================= Enumeration of k-chains of paths in the Tamari lattice, and valid chains ========================= """

# Check that a chain is valid for the nu-Tamari order of each of its element.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
def isValidChain(chain):
    for i in range(len(chain)-2):
        for j in range(i+1, len(chain)-1):
            if not compareNuTamari(chain[i], chain[j], chain[j+1]):
                return false
    return true

# WEAK CHAINS

@cached_function
# Return all (m,n)-Tamari weak chains of length k.
def TamariWeakChains(m,n,k):
    return weakChains(TamariLattice(m,n), k, element_constructor=tuple)

@cached_function
# Return all (m,n)-Tamari weak chains of length k which are valid.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
# We also require that the bottommost path appear as the first path.
def validTamariWeakChains(m,n,k):
    bp = bottomPath(m,n)
    validChains = []
    for chain in TamariWeakChains(m,n,k-1):
        extChain = (bp,) + chain
        if isValidChain(extChain):
            validChains.append(extChain)
    return validChains

# STRICT CHAINS

@cached_function
# Return all (m,n)-Tamari strict chains which are valid and between the given bottom path bp and top path tp.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
# By default, we require that the bottom path appear as the first path, and that the top path appears at the last path.
# This default settings can be changed with the optional arguments includeBottommost and includeTopmost.
def validTamariStrictChainsBetweenPaths(m, n, bp, tp, includeBottom=true, includeTop=true):
    tl = TamariLattice(m,n)
    interval = tl.interval(bp, tp)
    subTamariLattice = tl.subposet(interval[int(includeBottom):len(interval)-int(includeTop)])
    validStrictChains = []
    for chain in subTamariLattice.chains(element_constructor=tuple):
        chain = (bp,)*(int(includeBottom)) + chain + (tp,)*(int(includeTop))
        if isValidChain((bottomPath(m, n),)*(1-int(includeBottom))+chain):
            validStrictChains.append(chain)
    return validStrictChains

@cached_function
# Return all (m,n)-Tamari strict chains which are valid and start with the given path bp.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
def validTamariStrictChainsFromPath(m, n, bp, **kwargs):
    return validTamariStrictChainsBetweenPaths(m, n, bp, topPath(m, n), **kwargs)

@cached_function
# Return all (m,n)-Tamari strict chains which are valid and end with the given path tp.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
# We also require that the bottommost path appear as the first path.
def validTamariStrictChainsToPath(m, n, tp, **kwargs):
    return validTamariStrictChainsBetweenPaths(m, n, bottomPath(m, n), tp, **kwargs)

@cached_function
# Return all (m,n)-Tamari strict chains which are valid.
# This means that for all i < j < k, (Pj, Pk) is an interval in the Pi-Tamari lattice.
# By default, we require that the bottommost path appear as the first path, but not that the topmost path appears at the last path.
# This default settings can be changed with the optional arguments includeBottom and includeTop.
def validTamariStrictChains(m, n, includeBottom=true, includeTop=false):
    return validTamariStrictChainsBetweenPaths(m, n, bottomPath(m, n), topPath(m, n), includeBottom=includeBottom, includeTop=includeTop)

# MAXIMAL CHAINS

# Return the inclusion maximal chains among a given set of chains.
def maximalHopfChains(m, n):
    return Poset([validTamariStrictChains(m,n, includeBottom=true, includeTop=true), is_sublist]).maximal_elements()

""" ========================= Anklette and collar statistic ========================= """

# ANKLETTE STATISTIC
# Return the maximal Hopf distance between paths chain[-2] and chain[-1] with respect to chain.
# This means the length of the maximal (m,n)-Tamari chain between chain[-2] and chain[-1] which is also valid with respect to chain[:-2].
# This is called anklet statistic in the paper. 
def maximalHopfDistanceFirst(m,n,chain):
    if len(chain) == 1:
        return 0
    return max([len(superChain)-1 for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[:2], element_constructor=tuple) if isValidChain(superChain + chain[2:])])

# COLLAR STATISTIC
# Return the maximal Hopf distance between paths chain[-2] and chain[-1] with respect to chain.
# This means the length of the maximal (m,n)-Tamari chain between chain[-2] and chain[-1] which is also valid with respect to chain[:-2].
# This is called collar statistic in the paper. 
def maximalHopfDistanceLast(m,n,chain):
    if len(chain) == 1:
        return 0
    return max([len(superChain)-1 for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[-2:], element_constructor=tuple) if isValidChain(chain[:-2] + superChain)])

# DECOMPOSE THE ANKLETTE STATISTIC FOR THE DELTA CONJECTURE
# Return the level of the flip (row or column depending on direction).
def flipLevel(p, q, direction='C'):
    level = 0;
    for i in range(len(p)):
        if direction == 'C':
            if p[len(p)-1-i] != q[len(p)-1-i]:
                return level
            level += 1 - p[len(p)-1-i]
        if direction == 'R':
            if p[i] != q[i]:
                return level
            level += p[i]
    return -1

# Return the levels of a chain of flips (rows or columns depending on direction).
def flipLevels(m, n, chain, direction='C'):
    nlevels = m if direction == 'C' else n
    levels = [0]*nlevels
    for i in range(len(chain)-1):
        levels[flipLevel(chain[i], chain[i+1], direction=direction)] += 1
    return tuple(levels)

# Return the maximal Hopf distance by levels between paths chain[-2] and chain[-1] with respect to chain.
# This means the length of the maximal (m,n)-Tamari chain between chain[-2] and chain[-1] which is also valid with respect to chain[:-2].
# The levels are rows or columns depending on direction.
def maximalHopfDistanceLevels(m, n, chain, direction='C'):
    if len(chain) == 1:
        return ((0,)*m,)
    max = 0
    res = set([])
    for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[:2], element_constructor=tuple):
        if isValidChain(superChain + chain[2:]):
            if len(superChain)-1 > max:
                max = len(superChain)-1
                res = set([])
            if len(superChain)-1 == max:
                res.add(flipLevels(m, n, superChain, direction=direction))
    if len(res) > 1:
        print "I did not expect that..."
    return list(res)

# Return the height of an element x in a poset P, ie. the length of the heighest chain to this element.
def height(P, x):
    return P.subposet(P.principal_lower_set(x)).height()-1

# Compute the maximal Hopf distances in chain.
# For all i, we compute the maximal length of a superchain between chain[i] and chain[i+1] that is valid with the rest of chain.
def maximalHopfDistances(m,n,chain):
    if len(chain) == 1:
        return (0,)
    return tuple([max([len(superChain)-1 for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[i:i+2], element_constructor=tuple) if isValidChain(chain[:i] + superChain + chain[i+2:])]) for i in range(len(chain)-1)])

# Compute the maximal Hopf distances between the first two paths and the last two paths.
# This is just to gain time from the previous function
def maximalHopfDistancesFirstLast(m,n,chain):
    if len(chain) == 1:
        return (0,)
    return tuple([max([len(superChain)-1 for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[i:i+2], element_constructor=tuple) if isValidChain(chain[:i] + superChain + chain[i+2:])]) for i in [0, len(chain)-2]])

# Return all valid superchains of chain, with the tuple of lengths of the intervals between the elements of the initial chain.
def allValidSuperChainsWithLengthIntervals(m,n,chain):
    allValidSuperChains = []
    for chainsIntervals in cartesian_product(strictChainsIntervals(TamariLattice(m,n), chain)):
        lengthIntervals = [len(c)+1 for c in chainsIntervals]
        superChain = constructChain(chain, chainsIntervals)
        if isValidChain(superChain):
            allValidSuperChains.append((superChain, lengthIntervals))
    return allValidSuperChains

# Return the length and number of saturated valid chains containing the chain.
def recordMaximaLengthStrictSuperChains(m, n, chain):
    return recordMaxima(map(len, [superchain for superchain in strictSuperChains(TamariLattice(m,n), chain + (topPath(m,n),), element_constructor=tuple) if isValidChain(superchain)]))

""" ========================= Dinv of parking function and LLT polynomials ========================= """

# Transform a multipermutation to a parking function.
def pf2PF(path, pf):
    return ParkingFunction(labelling=Word(pf).standard_permutation().inverse(), area_sequence=path.to_area_sequence())

# Return the llt of a path, that is the sum of t^dinv(pf) * F_comp(pf).
# We use the function characteristic_quasisymmetric_function() already in sage
def llt1(path):
    return sum([pf2PF(path, pf).characteristic_quasisymmetric_function(R=R, q=R('q')) for pf in DyckWord(list(path)).list_parking_functions()])

"""
All what comes next is the same but avoids pasing through buildin sage functions.
It is just for me to check that I have the right definitions.
"""

# Return the diagonal levels of the north steps of a Dyck path.
def diagonalLevelsNorthSteps(path):
    diagonalLevels = []
    level = 0
    for x in path:
        if x == 1:
            diagonalLevels.append(level)
            level += 1
        if x == 0:
            level -= 1
    return diagonalLevels

# Return the positions of the north steps of path read by increasing levels and from south-west to north-east in each diagonal.
def levelReadingNorthSteps(path):
    return Word(diagonalLevelsNorthSteps(path)).standard_permutation().inverse()

# Return the dinv of a parking function pf over a Dyck path path.
# This is the number of diagonal inversions of the parking function on the Dyck path.
# A diagonal inversion is a pair of values i < j such that
#    - either these two values are on the same diagonal and i is south-west of j,
#    - or i is on the diagonal just below j and i is north-east of j.
def diagonalInversionsParkingFunction(path, pf):
    n = len(pf)
    # diagonalInversions = []
    numberDiagonalInversions = 0
    whichNorthSteps = Word(pf).standard_permutation() # which north step at position i in pf
    diagonalLevels = diagonalLevelsNorthSteps(path) # diagonal levels of the north steps of path
    diagonalLevelsReordered = [diagonalLevels[i-1] for i in whichNorthSteps]
    for i in range(n):
        for j in range(i+1,n):
            # easy case
            if diagonalLevelsReordered[i] == diagonalLevelsReordered[j] and pf[i] < pf[j]:
                # diagonalInversions.append([i,j])
                numberDiagonalInversions += 1
            # mobius case
            if diagonalLevelsReordered[i] == diagonalLevelsReordered[j]-1 and pf[i] > pf[j]:
                # diagonalInversions.append([i,j])
                numberDiagonalInversions += 1
    # return diagonalInversions
    return numberDiagonalInversions

# Return the permutation associated to a parking function as follows.
# We read the parking function by decreasing levels of the Dyck path, and in each level from north-east to south-west. 
# for instance for pf = [1, 1, 3, 3, 1, 3] we obtain [6, 4, 5, 3, 2, 1].
def permutationParkingFunction(pf):
    whichNorthStepsInverse = Word(pf).standard_permutation().inverse()
    return Permutation(reversed([whichNorthStepsInverse[i-1] for i in levelReadingNorthSteps(path)]))

def partitionParkingFunction(pf):
    return permutationParkingFunction(pf).inverse().descents_composition()

# return the llt of a path, that is the sum of the dinv of the parking functions over that path
def llt2(path):
    return sum([q^diagonalInversionsParkingFunction(path, pf) * F(partitionParkingFunction(pf)) for pf in DyckWord(list(path)).list_parking_functions()])

""" ========================= Check the enumeration ========================= """
    
@cached_function
# Refine the count according to the shape of the last path.
def refinedCountValidChains(m,n,k):
    validChains = validTamariWeakChains(m,n,k)
    res = dict({})
    for chain in validChains:
        upSteps = upStepsPartition(chain[k-1])
        if not res.has_key(upSteps):
            res[upSteps] = 0
        res[upSteps] += 1
    return res

# Return the symmetric polynomial expression.
def refinedCountValidChainsPolyn(m,n,k):
    return add([coeff * e(partition) for (partition, coeff) in  refinedCountValidChains(m,n,k).items()])

"""

sage: for n in range(2,7):
....:         print n, refinedCountValidChainsPolyn(3, n, 2) == e(Eval1(Phi[(3, n)], 2))

"""

""" ======================= r ======================= """

@cached_function
# Count the number of valid Tamari strict chains, according to the shape of the last path, and the length of the chain.
def refinedCountValidChains_r(m,n):
    validStrictChains = validTamariStrictChains(m,n)
    res = dict({})
    for chain in validStrictChains:
        key = (upStepsPartition(chain[-1]), len(chain))
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# Return the polynomial expression.
def refinedCountValidChainsPolyn_r(m,n):
    return add([coeff * binomial(r-1, length-1) * e(partition) for ((partition, length), coeff) in refinedCountValidChains_r(m,n).items()])

"""

sage: for m in range(2,7):
....:     for n in range(2,7):
....:         print m, n,  refinedCountValidChainsPolyn_r(m, n) == e(Eval1(Phi[(m, n)], r, {r}))
....:         
2 2 True
2 3 True
2 4 True
2 5 True
2 6 True
3 2 True
3 3 True
3 4 True
3 5 True
3 6 True
4 2 True
4 3 True
4 4 True
4 5 True
4 6 False
5 2 True
5 3 True
5 4 True
5 5 False
5 6 False
6 2 True
6 3 True
6 4 True
6 5 False
6 6 False

============================================== 55 ==============================================

sage: excess55 =  refinedCountValidChainsPolyn_r(5,5) - Eval1(Phi[(5,5)], r, {r})
sage: excess55
(1/5040*r^7-1/360*r^5+7/720*r^3-1/140*r)*e[3, 2] + (1/40320*r^8-1/10080*r^7-1/2880*r^6+1/720*r^5+7/5760*r^4-7/1440*r^3-1/1120*r^2+1/280*r)*e[4, 1] + (1/362880*r^9-1/12096*r^7+13/17280*r^5-41/18144*r^3+1/630*r)*e[5]
sage: expressSymFunctOnBinomialCoeffs(excess55)
[([3, 2], [0, 0, 0, 0, 1, 3, 3, 1]),
([4, 1], [0, 0, 0, 0, 0, 1, 3, 3, 1]),
([5], [0, 0, 0, 0, 0, 1, 4, 6, 4, 1])]


============================================== 46 ==============================================

sage: excess46 =  refinedCountValidChainsPolyn_r(4,6) - Eval1(Phi[(4,6)], r, {r})
sage: excess46
(1/24*r^4-1/4*r^3+11/24*r^2-1/4*r)*e[3, 3] + (1/120*r^5-1/12*r^4+7/24*r^3-5/12*r^2+1/5*r)*e[4, 2] + (1/720*r^6-1/80*r^5+5/144*r^4-1/48*r^3-13/360*r^2+1/30*r)*e[5, 1] + (1/5040*r^7-1/720*r^6+1/720*r^5+1/144*r^4-1/90*r^3-1/180*r^2+1/105*r)*e[6]
sage: expressSymFunctOnBinomialCoeffs(excess46)
[([3, 3], [0, 0, 0, 0, 1]),
([4, 2], [0, 0, 0, 0, 0, 1]),
([5, 1], [0, 0, 0, 0, 0, 1, 1]),
([6], [0, 0, 0, 0, 0, 1, 2, 1])]


============================================== 56 ==============================================

sage: excess56 =  refinedCountValidChainsPolyn_r(5,6) - Eval1(Phi[(5,6)], r, {r})
sage: excess56
(1/720*r^6-1/240*r^5-1/144*r^4+1/48*r^3+1/180*r^2-1/60*r)*e[3, 3] + (1/5040*r^7-1/720*r^6+1/720*r^5+1/144*r^4-1/90*r^3-1/180*r^2+1/105*r)*e[4, 2] + (1/40320*r^8-1/10080*r^7-1/2880*r^6+1/720*r^5+7/5760*r^4-7/1440*r^3-1/1120*r^2+1/280*r)*e[5, 1] + (1/362880*r^9-1/12096*r^7+13/17280*r^5-41/18144*r^3+1/630*r)*e[6]
sage: expressSymFunctOnBinomialCoeffs(excess56)
[([3, 3], [0, 0, 0, 0, 1, 2, 1]),
([6], [0, 0, 0, 0, 0, 1, 4, 6, 4, 1]),
([4, 2], [0, 0, 0, 0, 0, 1, 2, 1]),
([5, 1], [0, 0, 0, 0, 0, 1, 3, 3, 1])]


============================================== 65 ==============================================

sage: excess65 =  refinedCountValidChainsPolyn_r(6,5) - Eval1(Phi[(6,5)], r, {r})
sage: excess65
(1/5040*r^7-1/360*r^5+7/720*r^3-1/140*r)*e[3, 2] + (1/40320*r^8-1/10080*r^7-1/2880*r^6+1/720*r^5+7/5760*r^4-7/1440*r^3-1/1120*r^2+1/280*r)*e[4, 1] + (1/362880*r^9-1/12096*r^7+13/17280*r^5-41/18144*r^3+1/630*r)*e[5]
sage: expressSymFunctOnBinomialCoeffs(excess65)
[([3, 2], [0, 0, 0, 0, 1, 3, 3, 1]),
([4, 1], [0, 0, 0, 0, 0, 1, 3, 3, 1]),
([5], [0, 0, 0, 0, 0, 1, 4, 6, 4, 1])]


============================================== 66 ==============================================

sage: excess66 =  refinedCountValidChainsPolyn_r(6,6) - Eval1(Phi[(6,6)], r, {r})
sage: excess66
...
sage: expressSymFunctOnBinomialCoeffs(excess66)
...

"""

""" ======================= r q ======================= """

@cached_function
# Count the number of valid (m,n)-Tamari strict chains, according to the shape of the last path, the length of the chain, and the maximal Hopf distances between any two consecutive paths.
def refinedCountValidChains_r_q(m,n):
    validStrictChains = validTamariStrictChains(m,n)
    res = dict({})
    for chain in validStrictChains:
        key = (upStepsPartition(chain[-1]), len(chain), maximalHopfDistancesFirstLast(m,n,chain))
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# Return the symmetric polynomial expression, according to the shape of the last path and the maximal Hopf distance between the first two paths.
# This is the anklette statistic and it works (see below).
def refinedCountValidChainsPolyn_r_q_first(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) * e(partition) for ((partition, length, distances), coeff) in refinedCountValidChains_r_q(m,n).items()])

# Return the symmetric polynomial expression, according to the shape of the last path and the maximal Hopf distance between the last two paths.
# This is the collar statistic and it works (see below).
def refinedCountValidChainsPolyn_r_q_last(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * e(partition) for ((partition, length, distances), coeff) in refinedCountValidChains_r_q(m,n).items()])

"""

sage: for m in range(2,7):
....:     for n in range(2,7):
....:         Pmn = refinedCountValidChainsPolyn_r_q_first(m, n)
....:         Qmn = refinedCountValidChainsPolyn_r_q_last(m, n)
....:         Rmn = e(Eval1(Phi[(m, n)], r-1+q, {r}))
....:         print m, n, Pmn - Qmn, Pmn - Rmn, Qmn - Rmn
....:         
2 2 0 0 0
2 3 0 0 0
2 4 0 0 0
2 5 0 0 0
2 6 0 0 0
3 2 0 0 0
3 3 0 0 0
3 4 0 0 0
3 5 0 0 0
3 6 0 0 0
4 2 0 0 0
4 3 0 0 0
4 4 0 0 0
4 5 0 0 0
4 6
(-1/2*r^2*q^3+1/2*r^2*q^2+5/2*r*q^3-5/2*r*q^2-3*q^3+3*q^2)*e[4, 2] + (-1/6*r^3*q^3-1/2*r^2*q^4+1/6*r^3*q^2+2*r^2*q^3+5/2*r*q^4-3/2*r^2*q^2-41/6*r*q^3-3*q^4+13/3*r*q^2+7*q^3-4*q^2)*e[5, 1] + (-1/24*r^4*q^3-1/6*r^3*q^4-1/2*r^2*q^5+1/24*r^4*q^2+7/12*r^3*q^3+2*r^2*q^4+5/2*r*q^5-5/12*r^3*q^2-71/24*r^2*q^3-41/6*r*q^4-3*q^5+35/24*r^2*q^2+77/12*r*q^3+7*q^4-25/12*r*q^2-5*q^3+q^2)*e[6]
(1/24*r^4+1/6*r^3*q+1/2*r^2*q^2-5/12*r^3-3/2*r^2*q-5/2*r*q^2+35/24*r^2+13/3*r*q+3*q^2-25/12*r-4*q+1)*e[3, 3] + (1/120*r^5+1/24*r^4*q+1/6*r^3*q^2-1/8*r^4-7/12*r^3*q-3/2*r^2*q^2+17/24*r^3+71/24*r^2*q+13/3*r*q^2-15/8*r^2-77/12*r*q-4*q^2+137/60*r+5*q-1)*e[4, 2] + (1/720*r^6+1/120*r^5*q+1/24*r^4*q^2-1/48*r^5-1/8*r^4*q-5/12*r^3*q^2+17/144*r^4+17/24*r^3*q+35/24*r^2*q^2-5/16*r^3-15/8*r^2*q-25/12*r*q^2+137/360*r^2+137/60*r*q+q^2-1/6*r-q)*e[5, 1] + (1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2-1/360*r^6-1/48*r^5*q-1/12*r^4*q^2+1/72*r^5+17/144*r^4*q+7/24*r^3*q^2-1/36*r^4-5/16*r^3*q-5/12*r^2*q^2+7/720*r^3+137/360*r^2*q+1/5*r*q^2+11/360*r^2-1/6*r*q-1/42*r)*e[6]
(1/24*r^4+1/6*r^3*q+1/2*r^2*q^2-5/12*r^3-3/2*r^2*q-5/2*r*q^2+35/24*r^2+13/3*r*q+3*q^2-25/12*r-4*q+1)*e[3, 3] + (1/120*r^5+1/24*r^4*q+1/6*r^3*q^2+1/2*r^2*q^3-1/8*r^4-7/12*r^3*q-2*r^2*q^2-5/2*r*q^3+17/24*r^3+71/24*r^2*q+41/6*r*q^2+3*q^3-15/8*r^2-77/12*r*q-7*q^2+137/60*r+5*q-1)*e[4, 2] + (1/720*r^6+1/120*r^5*q+1/24*r^4*q^2+1/6*r^3*q^3+1/2*r^2*q^4-1/48*r^5-1/8*r^4*q-7/12*r^3*q^2-2*r^2*q^3-5/2*r*q^4+17/144*r^4+17/24*r^3*q+71/24*r^2*q^2+41/6*r*q^3+3*q^4-5/16*r^3-15/8*r^2*q-77/12*r*q^2-7*q^3+137/360*r^2+137/60*r*q+5*q^2-1/6*r-q)*e[5, 1] + (1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4+1/2*r^2*q^5-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-7/12*r^3*q^3-2*r^2*q^4-5/2*r*q^5+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+71/24*r^2*q^3+41/6*r*q^4+3*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-77/12*r*q^3-7*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+5*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[6]
5 2 0 0 0
5 3 0 0 0
5 4 0 0 0
5 5
(1/6*r^3*q^4-1/6*r^3*q^3-r^2*q^4+r^2*q^3+11/6*r*q^4-11/6*r*q^3-q^4+q^3)*e[3, 2] + (-1/6*r^3*q^4+1/6*r^3*q^3+r^2*q^4-r^2*q^3-11/6*r*q^4+11/6*r*q^3+q^4-q^3)*e[4, 1] + (-1/24*r^4*q^5+1/24*r^4*q^4+1/4*r^3*q^5-1/4*r^3*q^4-11/24*r^2*q^5+11/24*r^2*q^4+1/4*r*q^5-1/4*r*q^4)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/720*r^6-1/80*r^5*q-1/12*r^4*q^2-5/12*r^3*q^3-r^2*q^4+1/720*r^5+5/144*r^4*q+7/24*r^3*q^2+35/24*r^2*q^3+11/6*r*q^4+1/144*r^4-1/48*r^3*q-5/12*r^2*q^2-25/12*r*q^3-q^4-1/90*r^3-13/360*r^2*q+1/5*r*q^2+q^3-1/180*r^2+1/30*r*q+1/105*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/3360*r^7-1/360*r^6*q-1/48*r^5*q^2-1/8*r^4*q^3-5/12*r^3*q^4+1/960*r^6+1/72*r^5*q+17/144*r^4*q^2+17/24*r^3*q^3+35/24*r^2*q^4-1/36*r^4*q-5/16*r^3*q^2-15/8*r^2*q^3-25/12*r*q^4-11/1920*r^4+7/720*r^3*q+137/360*r^2*q^2+137/60*r*q^3+q^4+1/160*r^3+11/360*r^2*q-1/6*r*q^2-q^3+47/10080*r^2-1/42*r*q-1/168*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4-1/40320*r^8-1/3360*r^7*q-1/360*r^6*q^2-1/48*r^5*q^3-1/12*r^4*q^4+1/60480*r^7+1/960*r^6*q+1/72*r^5*q^2+17/144*r^4*q^3+7/24*r^3*q^4+1/2880*r^6-1/36*r^4*q^2-5/16*r^3*q^3-5/12*r^2*q^4-11/17280*r^5-11/1920*r^4*q+7/720*r^3*q^2+137/360*r^2*q^3+1/5*r*q^4-7/5760*r^4+1/160*r^3*q+11/360*r^2*q^2-1/6*r*q^3+59/22680*r^3+47/10080*r^2*q-1/42*r*q^2+1/1120*r^2-1/168*r*q-1/504*r)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/720*r^6-1/80*r^5*q-1/12*r^4*q^2-1/4*r^3*q^3+1/720*r^5+5/144*r^4*q+7/24*r^3*q^2+11/24*r^2*q^3+1/144*r^4-1/48*r^3*q-5/12*r^2*q^2-1/4*r*q^3-1/90*r^3-13/360*r^2*q+1/5*r*q^2-1/180*r^2+1/30*r*q+1/105*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/3360*r^7-1/360*r^6*q-1/48*r^5*q^2-1/8*r^4*q^3-1/4*r^3*q^4+1/960*r^6+1/72*r^5*q+17/144*r^4*q^2+13/24*r^3*q^3+11/24*r^2*q^4-1/36*r^4*q-5/16*r^3*q^2-7/8*r^2*q^3-1/4*r*q^4-11/1920*r^4+7/720*r^3*q+137/360*r^2*q^2+9/20*r*q^3+1/160*r^3+11/360*r^2*q-1/6*r*q^2+47/10080*r^2-1/42*r*q-1/168*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5-1/40320*r^8-1/3360*r^7*q-1/360*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5+1/60480*r^7+1/960*r^6*q+1/72*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5+1/2880*r^6-1/36*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5-11/17280*r^5-11/1920*r^4*q+7/720*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4-7/5760*r^4+1/160*r^3*q+11/360*r^2*q^2-1/6*r*q^3+59/22680*r^3+47/10080*r^2*q-1/42*r*q^2+1/1120*r^2-1/168*r*q-1/504*r)*e[5]

"""

# Return the symmetric polynomial expression, according to the shape of the last path and the maximal Hopf distance between the first two paths and between the last two paths. 
# This does not work: the anklette and collar statistics are not compatible together (see below).
def refinedCountValidChainsPolyn_r_q_firstLast(m,n):
    return add([coeff * (binomial(r-2, length-1) + (t^(distances[0]) + q^(distances[-1])) * binomial(r-3, length-2) + t^(distances[0]) * q^(distances[-1]) * binomial(r-3, length-3)) * e(partition) for ((partition, length, distances), coeff) in refinedCountValidChains_r_q(m,n).items()])

"""

sage: for i in range(2,7):
....:     for j in range(2,7):
....:         print i, j, refinedCountValidChainsPolyn_r_q_firstLast(i,j) - e(Eval1(Phi[(i,j)], r-2+q+t, {r}))
....:         
2 2 e[2]
2 3 e[3]
2 4 e[3, 1] + (r-2)*e[4]
2 5 e[4, 1] + (r-2)*e[5]
2 6 e[4, 2] + (r-2)*e[5, 1] + (1/2*r^2-3/2*r+1)*e[6]
3 2 e[2]
3 3 r*e[2, 1] + (1/2*r^2-1/2*r-2)*e[3]
3 4 e[2, 2] + (r-1)*e[3, 1] + (1/2*r^2-1/2*r-2)*e[4]
3 5 e[3, 1, 1] + (r-1)*e[3, 2] + (1/2*r^2+1/2*r-4)*e[4, 1] + (1/6*r^3+1/2*r^2-14/3*r+6)*e[5]
3 6 r*e[3, 2, 1] + (1/2*r^2-1/2*r-2)*e[3, 3] + (1/2*r^2-1/2*r-2)*e[4, 1, 1] + (q^2*t^2+1/6*r^3-q^2*t-q*t^2+1/2*r^2+q*t-11/3*r+5)*e[4, 2] + (q^3*t^2+1/24*r^4-q^3*t-q^2*t^2+7/12*r^3+q^2*t-61/24*r^2+11/12*r+3)*e[5, 1] + (q^4*t^2+1/120*r^5-q^4*t-q^3*t^2+1/6*r^4+q^3*t-17/24*r^3-2/3*r^2+26/5*r-5)*e[6]
4 2 e[1, 1] + (r-2)*e[2]
4 3 r*e[2, 1] + (1/2*r^2-1/2*r-2)*e[3]
4 4 (1/2*r^2+1/2*r)*e[2, 1, 1] + (1/6*r^3+1/2*r^2-8/3*r+1)*e[2, 2] + (q^2*t^3+1/24*r^4-q^2*t^2-q*t^3+7/12*r^3+q*t^2-37/24*r^2-1/12*r-1)*e[3, 1] + (q^3*t^3+1/120*r^5-q^3*t^2-q^2*t^3+1/6*r^4+q^2*t^2-13/24*r^3-2/3*r^2+31/30*r+3)*e[4]
4 5 r*e[2, 2, 1] + (1/2*r^2-1/2*r)*e[3, 1, 1] + (1/6*r^3+r^2-19/6*r-1)*e[3, 2] + (q^2*t^3+1/24*r^4-q^2*t^2-q*t^3+7/12*r^3+q*t^2-49/24*r^2+5/12*r+1)*e[4, 1] + (q^3*t^3+1/120*r^5-q^3*t^2-q^2*t^3+1/6*r^4+q^2*t^2-13/24*r^3-2/3*r^2+31/30*r+3)*e[5]
[...]
"""

""" ======================= r q t ======================= """

""" Careful: all this only works for m=n since I don't know how to deal with dinv in the rectangular case. """

@cached_function
# Count the number of valid (m,n)-Tamari strict chains, according to the shape of the last path, the length of the chain, the maximal Hopf distances between any two consecutive paths, and the dinv of the last path.
def refinedCountValidChains_r_q_t(m,n):
    validStrictChains = validTamariStrictChains(m,n)
    res = dict({})
    for chain in validStrictChains:
        key = (upStepsPartition(chain[-1]), len(chain), maximalHopfDistancesFirstLast(m,n,chain), DyckWord(list(chain[-1])).dinv())
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# Return the polynomial expression with distance between the first two paths.
def refinedCountValidChainsPolyn_r_q_1_first(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

# Return the polynomial expression with distance between the last two paths.
def refinedCountValidChainsPolyn_r_q_1_last(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

# Return the polynomial expression with dinv
def refinedCountValidChainsPolyn_r_1_q(m,n):
    return add([coeff * binomial(r-1, length-1) * q^dinv for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

"""
    
sage: for n in range(2,6):
....:     Pn = refinedCountValidChainsPolyn_r_q_1_first(n,n)
....:     Qn = refinedCountValidChainsPolyn_r_q_1_last(n,n)
....:     Rn = add([e(Eval1(Phi[(n,n)], r-1+q, {r})).coefficient(p) for p in Partitions(n)])
....:     print n, Pn - Qn, Pn - Rn, Qn - Rn
....:
2 0 0 0
3 0 0 0
4 0 0 0
5
-1/24*r^4*q^5 + 1/24*r^4*q^4 + 1/4*r^3*q^5 - 1/4*r^3*q^4 - 11/24*r^2*q^5 + 11/24*r^2*q^4 + 1/4*r*q^5 - 1/4*r*q^4
1/362880*r^9 + 1/40320*r^8*q + 1/5040*r^7*q^2 + 1/720*r^6*q^3 + 1/120*r^5*q^4 - 1/10080*r^7*q - 1/720*r^6*q^2 - 1/80*r^5*q^3 - 1/24*r^4*q^4 - 1/12096*r^7 - 1/2880*r^6*q + 1/720*r^5*q^2 + 5/144*r^4*q^3 + 1/24*r^3*q^4 + 1/720*r^5*q + 1/144*r^4*q^2 - 1/48*r^3*q^3 + 1/24*r^2*q^4 + 13/17280*r^5 + 7/5760*r^4*q - 1/90*r^3*q^2 - 13/360*r^2*q^3 - 1/20*r*q^4 - 7/1440*r^3*q - 1/180*r^2*q^2 + 1/30*r*q^3 - 41/18144*r^3 - 1/1120*r^2*q + 1/105*r*q^2 + 1/280*r*q + 1/630*r
1/362880*r^9 + 1/40320*r^8*q + 1/5040*r^7*q^2 + 1/720*r^6*q^3 + 1/120*r^5*q^4 + 1/24*r^4*q^5 - 1/10080*r^7*q - 1/720*r^6*q^2 - 1/80*r^5*q^3 - 1/12*r^4*q^4 - 1/4*r^3*q^5 - 1/12096*r^7 - 1/2880*r^6*q + 1/720*r^5*q^2 + 5/144*r^4*q^3 + 7/24*r^3*q^4 + 11/24*r^2*q^5 + 1/720*r^5*q + 1/144*r^4*q^2 - 1/48*r^3*q^3 - 5/12*r^2*q^4 - 1/4*r*q^5 + 13/17280*r^5 + 7/5760*r^4*q - 1/90*r^3*q^2 - 13/360*r^2*q^3 + 1/5*r*q^4 - 7/1440*r^3*q - 1/180*r^2*q^2 + 1/30*r*q^3 - 41/18144*r^3 - 1/1120*r^2*q + 1/105*r*q^2 + 1/280*r*q + 1/630*r

"""

# Return the polynomial expression with distance between the first two paths and between the last two paths.
def refinedCountValidChainsPolyn_r_q_t_first_last(m,n):
    return add([coeff * q^distances[0] * t^distances[-1] for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

""" check whether collar and anklet are equidistributed (yes) and whether the distribution is symmetric (no)

sage: for n in range(2,6):
....:     Pn = refinedCountValidChainsPolyn_r_q_t_first_last(n,n)
....:     Qn = Pn.subs({q:t, t:q})
....:     print n, Pn - Qn, Pn.subs({q:1}) - Pn.subs({t:1, q:t})
2 0 0
3 0 0
4 q^3*t^2 - q^2*t^3 - q^3*t + q*t^3 + q^2*t - q*t^2 0
"""

# Returns the polynomial expression distance between the first two paths with and dinv
def refinedCountValidChainsPolyn_r_q_t_first_flat(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) * t^dinv for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

# Returns the polynomial expression distance between the last two paths with and dinv
def refinedCountValidChainsPolyn_r_q_t_last_flat(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * t^dinv for ((partition, length, distances, dinv), coeff) in refinedCountValidChains_r_q_t(m,n).items()])

"""

sage: for n in range(2,6):
....:     Pn = refinedCountValidChainsPolyn_r_q_t_first_flat(n,n)
....:     Qn = Pn.subs({q:t, t:q})
....:     print n, Pn - Qn
....:
0
0
0
5 1/5040*r^7*q^2 + 1/720*r^6*q^3 + 1/120*r^5*q^4 + 1/720*r^6*q^2*t + 1/120*r^5*q^3*t + 1/24*r^4*q^4*t - 1/5040*r^7*t^2 - 1/720*r^6*q*t^2 - 1/720*r^6*t^3 - 1/120*r^5*q*t^3 - 1/120*r^5*t^4 - 1/24*r^4*q*t^4 - 1/5040*r^7*q - 1/360*r^6*q^2 - 1/48*r^5*q^3 - 1/12*r^4*q^4 + 1/5040*r^7*t - 1/80*r^5*q^2*t - 1/12*r^4*q^3*t - 1/4*r^3*q^4*t + 1/360*r^6*t^2 + 1/80*r^5*q*t^2 + 1/48*r^5*t^3 + 1/12*r^4*q*t^3 + 1/12*r^4*t^4 + 1/4*r^3*q*t^4 + 1/720*r^6*q + 1/72*r^5*q^2 + 17/144*r^4*q^3 + 7/24*r^3*q^4 - 1/720*r^6*t + 5/144*r^4*q^2*t + 7/24*r^3*q^3*t + 11/24*r^2*q^4*t - 1/72*r^5*t^2 - 5/144*r^4*q*t^2 - 17/144*r^4*t^3 - 7/24*r^3*q*t^3 - 7/24*r^3*t^4 - 11/24*r^2*q*t^4 - 1/720*r^5*q - 1/36*r^4*q^2 - 5/16*r^3*q^3 - 5/12*r^2*q^4 + 1/720*r^5*t - 1/48*r^3*q^2*t - 5/12*r^2*q^3*t - 1/4*r*q^4*t + 1/36*r^4*t^2 + 1/48*r^3*q*t^2 + 5/16*r^3*t^3 + 5/12*r^2*q*t^3 + 5/12*r^2*t^4 + 1/4*r*q*t^4 - 1/144*r^4*q + 7/720*r^3*q^2 + 137/360*r^2*q^3 + 1/5*r*q^4 + 1/144*r^4*t - 13/360*r^2*q^2*t + 1/5*r*q^3*t - 7/720*r^3*t^2 + 13/360*r^2*q*t^2 - 137/360*r^2*t^3 - 1/5*r*q*t^3 - 1/5*r*t^4 + 1/90*r^3*q + 11/360*r^2*q^2 - 1/6*r*q^3 - 1/90*r^3*t + 1/30*r*q^2*t - 11/360*r^2*t^2 - 1/30*r*q*t^2 + 1/6*r*t^3 + 1/180*r^2*q - 1/42*r*q^2 - 1/180*r^2*t + 1/42*r*t^2 - 1/105*r*q + 1/105*r*t

sage: for n in range(2,6):
....:     Pn = refinedCountValidChainsPolyn_r_q_t_first_flat(n,n)
....:     Qn = refinedCountValidChainsPolyn_r_q_t_last_flat(n,n)
....:     Rn = add([e(Eval1(Phi[(n,n)], r-2+q+t, {r})).coefficient(p) for p in Partitions(n)])
....:     print n, Pn - Qn, Pn - Rn, Qn - Rn
....:
2 0 0 0
3 0 0 0
4 0 0 0
5
-1/24*r^4*q^5 + 1/24*r^4*q^4 + 1/4*r^3*q^5 - 1/4*r^3*q^4 - 11/24*r^2*q^5 + 11/24*r^2*q^4 + 1/4*r*q^5 - 1/4*r*q^4
1/362880*r^9 + 1/40320*r^8*q + 1/5040*r^7*q^2 + 1/720*r^6*q^3 + 1/120*r^5*q^4 + 1/40320*r^8*t + 1/5040*r^7*q*t + 1/720*r^6*q^2*t + 1/120*r^5*q^3*t + 1/24*r^4*q^4*t - 1/40320*r^8 - 1/3360*r^7*q - 1/360*r^6*q^2 - 1/48*r^5*q^3 - 1/12*r^4*q^4 - 1/10080*r^7*t - 1/720*r^6*q*t - 1/80*r^5*q^2*t - 1/12*r^4*q^3*t - 1/4*r^3*q^4*t + 1/60480*r^7 + 1/960*r^6*q + 1/72*r^5*q^2 + 17/144*r^4*q^3 + 7/24*r^3*q^4 - 1/2880*r^6*t + 1/720*r^5*q*t + 5/144*r^4*q^2*t + 7/24*r^3*q^3*t + 11/24*r^2*q^4*t + 1/2880*r^6 - 1/36*r^4*q^2 - 5/16*r^3*q^3 - 5/12*r^2*q^4 + 1/720*r^5*t + 1/144*r^4*q*t - 1/48*r^3*q^2*t - 5/12*r^2*q^3*t - 1/4*r*q^4*t - 11/17280*r^5 - 11/1920*r^4*q + 7/720*r^3*q^2 + 137/360*r^2*q^3 + 1/5*r*q^4 + 7/5760*r^4*t - 1/90*r^3*q*t - 13/360*r^2*q^2*t + 1/5*r*q^3*t - 7/5760*r^4 + 1/160*r^3*q + 11/360*r^2*q^2 - 1/6*r*q^3 - 7/1440*r^3*t - 1/180*r^2*q*t + 1/30*r*q^2*t + 59/22680*r^3 + 47/10080*r^2*q - 1/42*r*q^2 - 1/1120*r^2*t + 1/105*r*q*t + 1/1120*r^2 - 1/168*r*q + 1/280*r*t - 1/504*r
1/362880*r^9 + 1/40320*r^8*q + 1/5040*r^7*q^2 + 1/720*r^6*q^3 + 1/120*r^5*q^4 + 1/24*r^4*q^5 + 1/40320*r^8*t + 1/5040*r^7*q*t + 1/720*r^6*q^2*t + 1/120*r^5*q^3*t + 1/24*r^4*q^4*t - 1/40320*r^8 - 1/3360*r^7*q - 1/360*r^6*q^2 - 1/48*r^5*q^3 - 1/8*r^4*q^4 - 1/4*r^3*q^5 - 1/10080*r^7*t - 1/720*r^6*q*t - 1/80*r^5*q^2*t - 1/12*r^4*q^3*t - 1/4*r^3*q^4*t + 1/60480*r^7 + 1/960*r^6*q + 1/72*r^5*q^2 + 17/144*r^4*q^3 + 13/24*r^3*q^4 + 11/24*r^2*q^5 - 1/2880*r^6*t + 1/720*r^5*q*t + 5/144*r^4*q^2*t + 7/24*r^3*q^3*t + 11/24*r^2*q^4*t + 1/2880*r^6 - 1/36*r^4*q^2 - 5/16*r^3*q^3 - 7/8*r^2*q^4 - 1/4*r*q^5 + 1/720*r^5*t + 1/144*r^4*q*t - 1/48*r^3*q^2*t - 5/12*r^2*q^3*t - 1/4*r*q^4*t - 11/17280*r^5 - 11/1920*r^4*q + 7/720*r^3*q^2 + 137/360*r^2*q^3 + 9/20*r*q^4 + 7/5760*r^4*t - 1/90*r^3*q*t - 13/360*r^2*q^2*t + 1/5*r*q^3*t - 7/5760*r^4 + 1/160*r^3*q + 11/360*r^2*q^2 - 1/6*r*q^3 - 7/1440*r^3*t - 1/180*r^2*q*t + 1/30*r*q^2*t + 59/22680*r^3 + 47/10080*r^2*q - 1/42*r*q^2 - 1/1120*r^2*t + 1/105*r*q*t + 1/1120*r^2 - 1/168*r*q + 1/280*r*t - 1/504*r

"""


""" ======================= r q and LLT ======================= """

@cached_function
# Count the number of valid (m,n)-Tamari strict chains, according to the last path, the length of the chain, the maximal Hopf distances between any two consecutive paths.
def refinedCountValidChains_r_q_LLT(m,n):
    validStrictChains = validTamariStrictChains(m,n)
    res = dict({})
    for chain in validStrictChains:
        key = (chain[-1], len(chain), maximalHopfDistancesFirstLast(m,n,chain))
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# Return the polynomial expression with the maximal Hopf distance between the first two paths.
def refinedCountValidChainsPolyn_r_q_LLT_first(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) * LLT[partitionPath(path)[:-1]] for ((path, length, distances), coeff) in refinedCountValidChains_r_q_LLT(m,n).items()])

# Return the polynomial expression with the maximal Hopf distance between the last two paths.
def refinedCountValidChainsPolyn_r_q_LLT_last(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * LLT[partitionPath(path)[:-1]] for ((path, length, distances), coeff) in refinedCountValidChains_r_q_LLT(m,n).items()])

# The next function is just used to print nice formulas for the paper.
# Return the polynomial expression with the maximal Hopf distance between the last two paths written in latex.
def refinedCountValidChainsPolyn_r_q_LLT_last_latex(m,n):
    return [(mu, latex(add(coeff * abstractBinomial(r-2,i) for (i,coeff) in enumerate(coeffs)))) for (mu,coeffs) in expressSymFunctOnBinomialCoeffs(e(refinedCountValidChainsPolyn_r_q_LLT_last(n,n)),-2)]


"""

sage: for n in range(2,6):
....:     Pn = refinedCountValidChainsPolyn_r_q_LLT_first(n,n)
....:     Qn = refinedCountValidChainsPolyn_r_q_LLT_last(n,n)
....:     Rn = Eval1(Phi[(n,n)], q+t+r-2, {r})
....:     print n, Pn - Qn, Pn - Rn, Qn - Rn
....:
2 0 0 0
3 0 0 0
4 0 0 0
5
(-1/24*r^4*q^5+1/24*r^4*q^4+1/4*r^3*q^5-1/4*r^3*q^4-11/24*r^2*q^5+11/24*r^2*q^4+1/4*r*q^5-1/4*r*q^4)*s[1, 1, 1, 1, 1] + (1/6*r^3*q^4-1/6*r^3*q^3-r^2*q^4+r^2*q^3+11/6*r*q^4-11/6*r*q^3-q^4+q^3)*s[2, 2, 1]
(1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/40320*r^8*t+1/5040*r^7*q*t+1/720*r^6*q^2*t+1/120*r^5*q^3*t+1/24*r^4*q^4*t-1/40320*r^8-1/3360*r^7*q-1/360*r^6*q^2-1/48*r^5*q^3-1/12*r^4*q^4-1/10080*r^7*t-1/720*r^6*q*t-1/80*r^5*q^2*t-1/12*r^4*q^3*t-1/4*r^3*q^4*t+1/60480*r^7+1/960*r^6*q+1/72*r^5*q^2+17/144*r^4*q^3+7/24*r^3*q^4-1/2880*r^6*t+1/720*r^5*q*t+5/144*r^4*q^2*t+7/24*r^3*q^3*t+11/24*r^2*q^4*t+1/2880*r^6-1/36*r^4*q^2-5/16*r^3*q^3-5/12*r^2*q^4+1/720*r^5*t+1/144*r^4*q*t-1/48*r^3*q^2*t-5/12*r^2*q^3*t-1/4*r*q^4*t-11/17280*r^5-11/1920*r^4*q+7/720*r^3*q^2+137/360*r^2*q^3+1/5*r*q^4+7/5760*r^4*t-1/90*r^3*q*t-13/360*r^2*q^2*t+1/5*r*q^3*t-7/5760*r^4+1/160*r^3*q+11/360*r^2*q^2-1/6*r*q^3-7/1440*r^3*t-1/180*r^2*q*t+1/30*r*q^2*t+59/22680*r^3+47/10080*r^2*q-1/42*r*q^2-1/1120*r^2*t+1/105*r*q*t+1/1120*r^2-1/168*r*q+1/280*r*t-1/504*r)*s[1, 1, 1, 1, 1] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/10080*r^7-1/720*r^6*q-1/80*r^5*q^2-1/12*r^4*q^3-1/4*r^3*q^4-1/2880*r^6+1/720*r^5*q+5/144*r^4*q^2+7/24*r^3*q^3+11/24*r^2*q^4+1/720*r^5+1/144*r^4*q-1/48*r^3*q^2-5/12*r^2*q^3-1/4*r*q^4+7/5760*r^4-1/90*r^3*q-13/360*r^2*q^2+1/5*r*q^3-7/1440*r^3-1/180*r^2*q+1/30*r*q^2-1/1120*r^2+1/105*r*q+1/280*r)*s[2, 1, 1, 1] + (1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/720*r^6-1/80*r^5*q-1/12*r^4*q^2-5/12*r^3*q^3-r^2*q^4+1/720*r^5+5/144*r^4*q+7/24*r^3*q^2+35/24*r^2*q^3+11/6*r*q^4+1/144*r^4-1/48*r^3*q-5/12*r^2*q^2-25/12*r*q^3-q^4-1/90*r^3-13/360*r^2*q+1/5*r*q^2+q^3-1/180*r^2+1/30*r*q+1/105*r)*s[2, 2, 1]
(1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5+1/40320*r^8*t+1/5040*r^7*q*t+1/720*r^6*q^2*t+1/120*r^5*q^3*t+1/24*r^4*q^4*t-1/40320*r^8-1/3360*r^7*q-1/360*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5-1/10080*r^7*t-1/720*r^6*q*t-1/80*r^5*q^2*t-1/12*r^4*q^3*t-1/4*r^3*q^4*t+1/60480*r^7+1/960*r^6*q+1/72*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5-1/2880*r^6*t+1/720*r^5*q*t+5/144*r^4*q^2*t+7/24*r^3*q^3*t+11/24*r^2*q^4*t+1/2880*r^6-1/36*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5+1/720*r^5*t+1/144*r^4*q*t-1/48*r^3*q^2*t-5/12*r^2*q^3*t-1/4*r*q^4*t-11/17280*r^5-11/1920*r^4*q+7/720*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4+7/5760*r^4*t-1/90*r^3*q*t-13/360*r^2*q^2*t+1/5*r*q^3*t-7/5760*r^4+1/160*r^3*q+11/360*r^2*q^2-1/6*r*q^3-7/1440*r^3*t-1/180*r^2*q*t+1/30*r*q^2*t+59/22680*r^3+47/10080*r^2*q-1/42*r*q^2-1/1120*r^2*t+1/105*r*q*t+1/1120*r^2-1/168*r*q+1/280*r*t-1/504*r)*s[1, 1, 1, 1, 1] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/10080*r^7-1/720*r^6*q-1/80*r^5*q^2-1/12*r^4*q^3-1/4*r^3*q^4-1/2880*r^6+1/720*r^5*q+5/144*r^4*q^2+7/24*r^3*q^3+11/24*r^2*q^4+1/720*r^5+1/144*r^4*q-1/48*r^3*q^2-5/12*r^2*q^3-1/4*r*q^4+7/5760*r^4-1/90*r^3*q-13/360*r^2*q^2+1/5*r*q^3-7/1440*r^3-1/180*r^2*q+1/30*r*q^2-1/1120*r^2+1/105*r*q+1/280*r)*s[2, 1, 1, 1] + (1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/720*r^6-1/80*r^5*q-1/12*r^4*q^2-1/4*r^3*q^3+1/720*r^5+5/144*r^4*q+7/24*r^3*q^2+11/24*r^2*q^3+1/144*r^4-1/48*r^3*q-5/12*r^2*q^2-1/4*r*q^3-1/90*r^3-13/360*r^2*q+1/5*r*q^2-1/180*r^2+1/30*r*q+1/105*r)*s[2, 2, 1]

"""

# remark: it gives the same function if we reverse the path!!!!
def checkWeirdSymmetry(m,n):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) * LLT[partitionPath(path)[:-1]] for ((path, length, distances), coeff) in refinedCountValidChains_r_q_LLT(m,n).items()]) == add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * LLT[partitionPath(reversePath(path))[:-1]] for ((path, length, distances), coeff) in refinedCountValidChains_r_q_LLT(m,n).items()])

"""

In the square case, we observe a weird symmetry.
This symmetry of course makes absolutly no sense on the rectangular case.
We are not sure whether or not to expect that this symmetry should hold for any square.

sage: for n in range(2,6):
....:     print n, checkWeirdSymmetry(n,n)
....: 
2 True
3 True
4 True
5 False

"""

""" ======================= Delta conjecture ========================= """

# We try to decompose the collar statistic according to the levels.

""" First r by r """

# count the number of valid (m,n)-Tamari strict chains, according to the shape of the last path, the length of the chain, and the maximal Hopf distance levels of last two paths
def refinedCountValidChainsDeltaConj_q(m, n, r, direction='C'):
    (direction, reverseDirection) = ('C','R') if direction == 'C' else ('R','C')
    validWeakChains = validTamariWeakChains(m, n, r)
    upStepsHDCValidWeakChainsDict = dict({})
    for chain in validWeakChains:
        key = (upStepsPartition(chain[-1]), tuple(maximalHopfDistanceLevels(m, n, chain, direction=reverseDirection)), tuple(maximalHopfDistanceLevels(m, n, chain, direction=direction)))
        # print t, key
        if not upStepsHDCValidWeakChainsDict.has_key(key):
            upStepsHDCValidWeakChainsDict[key] = 0
        upStepsHDCValidWeakChainsDict[key] += 1
    return upStepsHDCValidWeakChainsDict

# returns the polynomial expression
def refinedCountValidChainsPolynDeltaConj_q(m, n, k, r, direction='C'):
    upStepsHDCValidWeakChainsDict = refinedCountValidChainsDeltaConj(m, n, r, direction=direction)
    res = 0
    for ((partition, distanceRows, distanceCols), coeff) in upStepsHDCValidWeakChainsDict.items():
        if len(distanceCols) == 1:
            distanceCols = distanceCols[0]
            nondescents = [i for i in range(1, len(distanceCols)) if distanceCols[i] > distanceCols[i-1]]
            # print partition, distanceRows, distanceCols, nondescents
            for s in Subsets(nondescents, k):
                # rem is a subset of size n-1-k containing all descents
                rem = set(range(1,len(distanceCols))).difference(s)
                # print rem, prod([q^distanceCols[n-i] for i in rem]) * e(partition)
                res += coeff * prod([q^distanceCols[i] for i in rem]) * e(partition)
    return res

def checkDeltaConjecture_q(m, n, k, r, direction='C'):
    A = refinedCountValidChainsPolynDeltaConj_q(m, n, k, r, direction=direction)
    B = e(Eval1(Skew1(e[k], Phi[(m,n)]), (r-1+q)*s[0]))
    return A == B

"""
    
sage: for n in range(2,7):
....:     for k in range(n):
....:         print n, k, checkDeltaConjecture(n, n, k, 2)
....: 
2 0 True
2 1 True
2 2 True
2 3 True
2 4 True
3 0 True
3 1 True
3 2 True
3 3 True
3 4 True
4 0 True
4 1 True
4 2 True
4 3 True
4 4 True
5 0 True
[...]

sage: for i in range(2,7):
....:     for k in range(0,6):
....:         print i, k, checkDeltaConjecture(i, i, k, 3)
....: 
2 0 True
2 1 True
2 2 True
2 3 True
2 4 True
3 0 True
3 1 True
3 2 True
3 3 True
3 4 True
4 0 True
4 1 True
4 2 True
4 3 True
4 4 True
5 0 I did not expect that...
I did not expect that...
True
5 1 I did not expect that...
I did not expect that...
False
5 2 I did not expect that...
I did not expect that...
False
5 3 I did not expect that...
I did not expect that...
True
5 4 I did not expect that...
I did not expect that...
True
[...]

sage: sage: for i in range(2,7):
....: ....:     for k in range(0,5):
....: ....:         print i, k, checkDeltaConjecture(i, i, k, 4)
....: 
2 0 True
2 1 True
2 2 True
2 3 True
2 4 True
3 0 True
3 1 True
3 2 True
3 3 True
3 4 True
4 0 True
4 1 True
4 2 True
4 3 True
4 4 True
[...]

"""

@cached_function
# count the number of valid (m,n)-Tamari strict chains, according to the shape of the last path, the length of the chain, and the maximal Hopf distance levels of last two paths
def refinedCountValidChainsDeltaConj(m, n, direction='C'):
    (direction, reverseDirection) = ('C','R') if direction == 'C' else ('R','C')
    validStrictChains = validTamariStrictChains(m, n)
    res = dict({})
    for chain in validStrictChains:
        key = (chain[-1], len(chain), tuple(maximalHopfDistanceLevels(m, n, chain, direction=reverseDirection)), tuple(maximalHopfDistanceLevels(m, n, chain, direction=direction)))
        # print t, key
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# returns the polynomial expression
def refinedCountValidChainsPolynDeltaConj_r_q(m, n, k, direction='C'):
    topPathLengthHDCValidWeakChainsDict = refinedCountValidChainsDeltaConj(m, n, direction=direction)
    res = 0
    for ((path, length, distanceRows, distanceCols), coeff) in topPathLengthHDCValidWeakChainsDict.items():
        if len(distanceCols) == 1:
            distanceCols = distanceCols[0]
            nondescents = [i for i in range(1, len(distanceCols)) if distanceCols[i] > distanceCols[i-1]]
            for s in Subsets(nondescents, k):
                # rem is a subset of size n-1-k containing all descents
                rem = set(range(1,len(distanceCols))).difference(s)
                res += coeff * (int(k==0) * binomial(r-2, length-1) + prod([q^distanceCols[i] for i in rem]) * binomial(r-2, length-2)) * e(upStepsPartition(path))
    return res

def checkDeltaConjecture_r_q(m, n, k, direction='C'):
    A = refinedCountValidChainsPolynDeltaConj_r_q(m, n, k, direction=direction)
    B = e(Eval1(Skew1(e[k], Phi[(m,n)]), r-1+q, {r}))
    return A == B

"""

sage: for n in range(2,7):
....:     for k in range(n):
....:         print n, k, checkDeltaConjecture_r_q(n, n, k)
....: 
2 0 True
2 1 True
2 2 True
2 3 True
2 4 True
3 0 True
3 1 True
3 2 True
3 3 True
3 4 True
4 0 True
4 1 True
4 2 True
4 3 True
4 4 True
5 0 I did not expect that...
[...]

Same for rectangular Delta conjecture:

sage: for m in range(2,7):
....:     for n in range(2,7):
....:         for k in range(n):
....:             print m, n, k, checkDeltaConjecture_r_q(m, n, k)
....: 
2 2 0 True
2 2 1 True
2 3 0 True
2 3 1 True
2 3 2 True
2 4 0 True
2 4 1 False
2 4 2 True
2 4 3 True
2 5 0 True
2 5 1 False
2 5 2 True
2 5 3 True
2 5 4 True
2 6 0 True
2 6 1 False
2 6 2 True
2 6 3 True
2 6 4 True
2 6 5 True
3 2 0 True
3 2 1 True
3 3 0 True
3 3 1 True
3 3 2 True
3 4 0 True
3 4 1 True
3 4 2 True
3 4 3 True
3 5 0 True
3 5 1 False
3 5 2 False
3 5 3 True
3 5 4 True
3 6 0 True
3 6 1 False
3 6 2 False
3 6 3 True
3 6 4 True
3 6 5 True
4 2 0 True
4 2 1 True
4 3 0 True
4 3 1 True
4 3 2 True
4 4 0 True
4 4 1 True
4 4 2 True
4 4 3 True
4 5 0 True
4 5 1 True
4 5 2 True
4 5 3 True
4 5 4 True
4 6 0 I did not expect that...
I did not expect that...
I did not expect that...
I did not expect that...
I did not expect that...
I did not expect that...
I did not expect that...
I did not expect that...
False
4 6 1 False
4 6 2 False
4 6 3 False
4 6 4 True
4 6 5 True
5 2 0 True
5 2 1 True
5 3 0 True
5 3 1 True
5 3 2 True
5 4 0 True
5 4 1 True
5 4 2 True
5 4 3 True

"""

# returns the polynomial expression
def refinedCountValidChainsPolynDeltaConj_r_q_LLT(m, n, k, direction='C'):
    topPathLengthHDCValidWeakChainsDict = refinedCountValidChainsDeltaConj(m, n, direction=direction)
    res = 0
    for ((path, length, distanceRows, distanceCols), coeff) in topPathLengthHDCValidWeakChainsDict.items():
        if len(distanceCols) == 1:
            distanceCols = distanceCols[0]
            nondescents = [i for i in range(1, len(distanceCols)) if distanceCols[i] > distanceCols[i-1]]
            for s in Subsets(nondescents, k):
                # rem is a subset of size n-1-k containing all descents
                rem = set(range(1,len(distanceCols))).difference(s)
                res += coeff * (int(k==0) * binomial(r-2, length-1) + prod([q^distanceCols[i] for i in rem]) * binomial(r-2, length-2)) * LLT[partitionPath(path)[:-1]]
    return res

def checkDeltaConjecture_r_q_LLT(m, n, k, direction='C'):
    A = refinedCountValidChainsPolynDeltaConj_r_q_LLT(m, n, k, direction=direction)
    B = e(Eval1(Skew1(e[k], Phi[(m,n)]), q+t+r-2, {r}))
    return A == B

"""
    
sage: for n in range(2,7):
....:     for k in range(n):
....:         print n, k, checkDeltaConjecture_r_q_LLT(n, n, k)
....: 
2 0 True
2 1 True
3 0 True
3 1 True
3 2 True
4 0 True
4 1 True
4 2 True
4 3 True
5 0 I did not expect that...
[...]

"""

""" ======================= The Cesar-Francois conjecture ========================= """

# compute the number of valid Hopf chains below all paths with a given up steps partition
def numberValidTamariStrictChainsToPartition(m, n, partition):
    res = 0
    for path in TamariLattice(m,n):
        if upStepsPartition(path) == partition:
            res += len(validTamariStrictChainsToPath(m, n, path))
    return res

# compute the sum of the coefficients in the M basis
def numberQSymMonomial(m, n, partition):
    res = 0
    for partition2, coeff in Sym.monomial()(Scalar2(Pleth1(Phi[(m,n)], s[1]+1), f(partition))):
        res += coeff * sum(QSym(Sym.monomial()(partition2)).coefficients())
    return res

# check that both coincide for all paths
def testNumbersMonomials(m, n):
    return all([numberValidTamariStrictChainsToPartition(m, n, partition) == numberQSymMonomial(m, n, partition) for partition in Partitions(n)])

def testNumbersMonomialsPrint(m, n):
    for partition in Partitions(n):
        print partition, numberValidTamariStrictChainsToPartition(m, n, partition) - numberQSymMonomial(m, n, partition)

"""

sage: testNumbersMonomialsPrint(3,3)
[3] 0
[2, 1] 0
[1, 1, 1] 0

sage: testNumbersMonomialsPrint(4,4)
[4] 0
[3, 1] 0
[2, 2] 0
[2, 1, 1] 0
[1, 1, 1, 1] 0

sage: testNumbersMonomialsPrint(3,4)
[4] 0
[3, 1] 0
[2, 2] 0
[2, 1, 1] 0
[1, 1, 1, 1] 0

sage: testNumbersMonomialsPrint(4,3)
[3] 0
[2, 1] 0
[1, 1, 1] 0

sage: testNumbersMonomialsPrint(3,6)
[6] 0
[5, 1] 0
[4, 2] 0
[4, 1, 1] 0
[3, 3] 0
[3, 2, 1] 0
[3, 1, 1, 1] 0
[2, 2, 2] 0
[2, 2, 1, 1] 0
[2, 1, 1, 1, 1] 0
[1, 1, 1, 1, 1, 1] 0

sage: for m in range(2,7):
....:     print m, testNumbersMonomials(m,3)
....:     
2 True
3 True
4 True
5 True
6 True

sage: for m in range(2,10):
....:     print m, testNumbersMonomials(m,4)
....:     
2 True
3 True
4 True
5 True
6 True
7 True
8 False
9 False

sage: for n in range(2,7):
....:     print n, testNumbersMonomials(4,n)
....:
2 True
3 True
4 True
5 True
6 False

"""

# Some data to help Francois 
def HopfChainsBelowPathCollarAnklette(m, n, path):
    resFirst = 0
    # resLast = 0
    for chain in validTamariStrictChainsToPath(m,n,path):
        resFirst += q^maximalHopfDistanceFirst(m,n,chain)
        #resLast += q^maximalHopfDistanceLast(m,n,chain)
    # if resFirst != resLast:
        #print "!!!! anklette != collar !!!!"
    return resFirst

# Some pictures to help François and Cesar
def TamariLatticeWithNumberHopfChains(m, n):
    tl = TamariLattice(m,n)
    dictPaths = {path:i for (i,path) in enumerate(tl)}
    return tl.relabel(lambda path: (dictPaths[path], HopfChainsBelowPathCollarAnklette(m,n,path).subs(q=1)))

def helpFrancois(n):
    for path in DyckWords(n-1):
        extPath = Word([1]+list(path)+[0])
        polExtPath =  HopfChainsBelowPathCollarAnklette(n, n, extPath)
        print partitionPath(extPath), "---", polExtPath, "---", polExtPath.subs(q=1)

# This was a conjecture of Francois...
def quasiSymCesarFrancois(m, n, path):
    return add([M(list(maximalHopfDistances(m, n, chain))) for chain in validTamariStrictChainsToPath(m,n,path)])

@cached_function
# Return a dictionnary dictChains such that dictChains[(mp,tp))] is the number of strict Hopf chains containing mp and finishing with tp. 
def dictMasterFormula(m, n):
    dictChains = defaultdict(lambda: 0)
    for chain in validTamariStrictChains(m, n):
        for path in chain:
            dictChains[(partitionPath(reversePath(path)), partitionPath(reversePath(chain[-1])))] += 1
    return dictChains

# Return two dictionnaries dictChainsMP and dictChainsTP such that dictChainsMP[mp][tp] and dictChainsTP[tp][mp] are both the number of strict Hopf chains containing mp and finishing with tp. 
def dictsMasterFormula(m, n):
    dictChainsMP = defaultdict(lambda: defaultdict(lambda: 0))
    dictChainsTP = defaultdict(lambda: defaultdict(lambda: 0))
    for (mp, tp), count in dictMasterFormula(m,n).items():
        dictChainsMP[mp][tp] = count
        dictChainsTP[tp][mp] = count
    return dictChainsMP, dictChainsTP

# nice print of the previous pair of dictionnaries
def printDictsMasterFormula(m, n):
    dictChainsMP, dictChainsTP = dictsMasterFormula(m, n)
    for path in dictChainsMP.keys():
        print path
        print add(dictChainsMP[path].values())
        print dictChainsMP[path].items()
        print add(dictChainsTP[path].values())
        print dictChainsTP[path].items()
        print

# returns a dictionnary dictChains such that dictChains[tp] is the number of strict Hopf chains starting with the staircase path and finishing with tp. 
@cached_function
def dictTP(m, n):
    dictChains = defaultdict(lambda: 0)
    for chain in validTamariStrictChains(m, n, includeBottom=true, includeTop=false):
        dictChains[partitionPath(reversePath(chain[-1]))] += 1
    return dictChains

# Returns a dictionnary dictChains such that dictChains[(mp,tp))] is the number of strict Hopf chains starting with mp and finishing with tp. 
@cached_function
def dictMPTP(m, n):
    dictChains = defaultdict(lambda: 0)
    for chain in validTamariStrictChains(m, n, includeBottom=false, includeTop=false):
        if len(chain) > 0:
            dictChains[(partitionPath(reversePath(chain[0])), partitionPath(reversePath(chain[-1])))] += 1
    return dictChains

# Return two dictionnaries dictChainsMP and dictChainsTP such that dictChainsMP[mp][tp] and dictChainsTP[tp][mp] are both the number of strict Hopf chains starting with mp and finishing with tp.
@cached_function
def dictsMPTP(m, n):
    dictChainsMP = defaultdict(lambda: defaultdict(lambda: 0))
    dictChainsTP = defaultdict(lambda: defaultdict(lambda: 0))
    for (mp, tp), count in dictMPTP(m,n).items():
        dictChainsMP[mp][tp] = count * dictTP(m, n)[mp]
        dictChainsTP[tp][mp] = count * dictTP(m, n)[mp]
    return dictChainsMP, dictChainsTP

# nice print of the previous pair of dictionnaries
def printDictsMPTP(m, n):
    dictChainsMP, dictChainsTP = dictsMPTP(m, n)
    for path in dictChainsMP.keys():
        print path
        print add(dictChainsMP[path].values())
        print dictChainsMP[path].items()
        print add(dictChainsTP[path].values())
        print dictChainsTP[path].items()
        print

"""

sage: for path in DyckWords(4):
....:     qsCF = quasiSymCesarFrancois(4,4,Word(path))
....:     print path
....:     if qsCF.is_symmetric():
....:         print s(qsCF.to_symmetric_function())
....:     else:
....:         print F(qsCF)
....:     
()()()()
s[]
()()(())
s[1]
()(())()
s[1]
()(()())
s[2]
()((()))
s[1, 1] + s[3]
(())()()
s[1]
(())(())
s[1, 1] + s[2]
(()())()
s[2]
(()()())
s[3]
(()(()))
s[2, 1] + s[4]
((()))()
s[1, 1] + s[3]
((())())
s[2, 1] + s[4]
((()()))
F[1, 3] + F[2, 1] + 2*F[2, 1, 1] + F[2, 1, 1, 1] - F[2, 1, 2] - F[3, 1, 1] + F[3, 2] + F[5]
(((())))
F[1, 1, 1] + F[1, 3] + F[1, 4] + F[2, 1, 1, 2] + 2*F[2, 1, 2] - F[2, 1, 3] + 2*F[2, 2] + F[3, 1] - F[3, 1, 2] + F[3, 3] + F[4, 1] + F[6]

"""

"""
sage: for path in DyckWords(4):
....:     print path, quasiSymCesarFrancois(4,4,Word(path))

"""

"""
sage: path = Word([1,1,1,0,1,0,0,0])                           
sage: allChainsToPath = validTamariStrictChainsToPath(4,4,path)
sage: for chain in allChainsToPath:
....:     print
....:     print
....:     print chain
....:     lengthsSuperChains = [l for (x,l) in allValidSuperChainsWithLengthIntervals(4,4,chain)]
....:     maxTotalLength = max(add(l) for l in lengthsSuperChains)
....:     print [l for l in lengthsSuperChains if add(l) == maxTotalLength]
....:     print maximalHopfDistances(4,4,chain)
....:     


(word: 10101010, word: 11101000)
[[5]]
(5,)


(word: 10101010, word: 10110010, word: 11101000)
[[1, 3]]
(1, 3)


(word: 10101010, word: 10110010, word: 11010010, word: 11101000)
[[1, 1, 1]]
(1, 1, 1)


(word: 10101010, word: 10110010, word: 10111000, word: 11101000)
[[1, 1, 2]]
(1, 1, 2)


(word: 10101010, word: 10110010, word: 10111000, word: 11011000, word: 11101000)
[[1, 1, 1, 1]]
(1, 1, 1, 1)


(word: 10101010, word: 10110010, word: 11011000, word: 11101000)
[[1, 2, 1]]
(1, 2, 1)


(word: 10101010, word: 11010010, word: 11101000)
[[2, 1]]
(2, 1)


(word: 10101010, word: 10101100, word: 11101000)
[[1, 4]]
(1, 4)


(word: 10101010, word: 10101100, word: 10110100, word: 11101000)
[[1, 1, 3]]
(1, 1, 3)


(word: 10101010, word: 10101100, word: 10110100, word: 10111000, word: 11101000)
[[1, 1, 1, 2]]
(1, 1, 1, 2)


(word: 10101010, word: 10101100, word: 10110100, word: 10111000, word: 11011000, word: 11101000)
[[1, 1, 1, 1, 1]]
(1, 1, 1, 1, 1)


(word: 10101010, word: 10101100, word: 10110100, word: 11010100, word: 11101000)
[[1, 1, 1, 1]]
(1, 1, 1, 1)


(word: 10101010, word: 10101100, word: 10110100, word: 11011000, word: 11101000)
[[1, 1, 2, 1]]
(1, 1, 2, 1)


(word: 10101010, word: 10101100, word: 10111000, word: 11101000)
[[1, 2, 2]]
(1, 2, 2)


(word: 10101010, word: 10101100, word: 10111000, word: 11011000, word: 11101000)
[[1, 2, 1, 1]]
(1, 2, 1, 1)


(word: 10101010, word: 10101100, word: 11010100, word: 11101000)
[[1, 1, 2], [1, 2, 1]]
(1, 2, 2)


(word: 10101010, word: 10101100, word: 11010100, word: 11011000, word: 11101000)
[[1, 1, 1, 1]]
(1, 1, 1, 1)


(word: 10101010, word: 10101100, word: 11011000, word: 11101000)
[[1, 3, 1]]
(1, 3, 1)


(word: 10101010, word: 10110100, word: 11101000)
[[2, 3]]
(2, 3)


(word: 10101010, word: 10110100, word: 10111000, word: 11101000)
[[2, 1, 2]]
(2, 1, 2)


(word: 10101010, word: 10110100, word: 10111000, word: 11011000, word: 11101000)
[[2, 1, 1, 1]]
(2, 1, 1, 1)


(word: 10101010, word: 10110100, word: 11010100, word: 11101000)
[[2, 1, 1]]
(2, 1, 1)


(word: 10101010, word: 10110100, word: 11011000, word: 11101000)
[[2, 2, 1]]
(2, 2, 1)


(word: 10101010, word: 10111000, word: 11101000)
[[3, 2]]
(3, 2)


(word: 10101010, word: 10111000, word: 11011000, word: 11101000)
[[3, 1, 1]]
(3, 1, 1)


(word: 10101010, word: 11010100, word: 11101000)
[[2, 2], [3, 1]]
(3, 2)


(word: 10101010, word: 11010100, word: 11011000, word: 11101000)
[[2, 1, 1]]
(2, 1, 1)


(word: 10101010, word: 11011000, word: 11101000)
[[4, 1]]
(4, 1)

"""

""" ======================= LOOKING FOR THE KILLER ========================= """

"""
    Now we try to understand the killers:
    - compute the chains that are killed by a given killer
    - find candidates for the killers
"""

def candidateKillers(m, n, goal, mu):
    killers = []
    bp = bottomPath(m,n) # the bottom path
    for tp in upStepDictionnary(m,n)[mu]: # the possible top paths
        belowtp = TamariLattice(m,n).subposet(TamariLattice(m,n).principal_lower_set(tp))
        for candidate in [c for c in belowtp.chains() if len(c) == 4 and c[0] == bp and c[-1] == tp]:
            contribution = sum(x^len(c) for c in strictSuperChains(belowtp, candidate, element_constructor=tuple) if isValidChain(c))
            if contribution == goal:
                killers.append(candidate)
    return killers

"""

============================================== 55 ==============================================

sage: candidates55 = candidateKillers(5, 5, add(x^len(c) for c in Subsets(range(3))) * x^4, Partition([3,2]))
sage: candidates55
[[word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000],
[word: 1010101010, word: 1010110100, word: 1011100100, word: 1101110000],
[word: 1010101010, word: 1011010100, word: 1011100100, word: 1101110000],
[word: 1010101010, word: 1011010100, word: 1101100100, word: 1110110000],
[word: 1010101010, word: 1101010100, word: 1101101000, word: 1110110000],
[word: 1010101010, word: 1010110100, word: 1011011000, word: 1110011000],
[word: 1010101010, word: 1010110100, word: 1101011000, word: 1110011000]]
sage: candidateKillers(5, 5, sum(x^len(c) for c in Subsets(range(3))) * x^5, Partition([4,1]))
[]


============================================== 46 ==============================================

sage: candidates46 = candidateKillers(4, 6, x^4, Partition([3,3]))
[(word: 1101011010, word: 1101101010, word: 1110110100, word: 1110111000)]


============================================== 56 ==============================================

sage: candidates56 = candidateKillers(5, 6, sum(x^len(c) for c in Subsets(range(2))) * x^4, Partition([3,3]))
[(word: 11010101010, word: 11011010010, word: 11011101000, word: 11101110000),
(word: 11010101010, word: 11101010100, word: 11101011000, word: 11101110000),
(word: 11010101010, word: 11101010100, word: 11101101000, word: 11101110000),
(word: 11010101010, word: 11010110010, word: 11011100010, word: 11101110000),
(word: 11010101010, word: 11011010010, word: 11011100010, word: 11101110000),
(word: 11010101010, word: 11010110010, word: 11101011000, word: 11101110000),
(word: 11010101010, word: 11011010010, word: 11011101000, word: 11101110000),
(word: 11010101010, word: 11011001010, word: 11011100100, word: 11101110000),
(word: 11010101010, word: 11011001010, word: 11011100100, word: 11101110000)]


============================================== 65 ==============================================

sage: candidates65 = candidateKillers(6, 5, sum(x^len(c) for c in Subsets(range(3))) * x^4, Partition([3,2]))
[...] (240 candidates)

"""

"""
Once we have candidates for killers, we can try to use fingerprints (such as the formula with the collar statistic, or with the anklet statistic, or with the dinv, or...) to choose among the killers.
"""

# Compute the maximal Hopf distances between the first two paths and the last two paths with respect to chain and a set of forbidden chains.
# This means the lengths of the maximal Tamari chain between first two paths and the last two paths which is also valid with respect to the rest of the chain and a set of forbidden chains
def maximalHopfDistancesFirstLastForbidden(m, n, chain, forbiddenChains):
    if len(chain) == 1:
        return (0,)
    return tuple([max([len(superChain)-1 for superChain in strictBoundedSuperChains(TamariLattice(m,n), chain[i:i+2], element_constructor=tuple) if isValidChain(chain[:i] + superChain + chain[i+2:]) and not chain[:i] + superChain + chain[i+2:] in forbiddenChains]) for i in [0, len(chain)-2]])

# Count the number of valid Tamari strict chains, according to the last path, the shape of the last path, the length of the chain, and the maximal Hopf distances between the first two paths and the last two paths with respect to a set of forbidden chains.
def refinedCountValidChainsFirstLastForbidden(m, n, forbiddenChains):
    validStrictChains = set(validTamariStrictChains(m,n)).difference(forbiddenChains)
    res = dict({})
    for chain in validStrictChains:
        key = (chain[-1], len(chain), maximalHopfDistancesFirstLastForbidden(m, n, chain, forbiddenChains))
        if not res.has_key(key):
            res[key] = 0
        res[key] += 1
    return res

# Return the polynomial expression.
def refinedCountValidChainsPolynForbidden_r(m, n, forbiddenChains):
    return add([coeff * binomial(r-1, length-1) * e(upStepsPartition(path)) for ((path, length, distances), coeff) in refinedCountValidChainsFirstLastForbidden(m, n, forbiddenChains).items()])

# choose among the candidates using the statistics q^collar and t^dinv
def chooseKiller1(m, n, candidates):
    goal = e(Eval1(Phi[(m,n)], r, {r}))
    keepCandidates = []
    for candidate in candidates:
        print
        print "---"
        print
        print candidate
        stats = refinedCountValidChainsPolynForbidden_r(m, n, strictSuperChains(TamariLattice(m,n), candidate, element_constructor=tuple))
        print stats - goal
        if stats == goal:
            keepCandidates.append(candidate)
            print "this is the killer"
    return keepCandidates

"""

============================================== 55 ==============================================

sage: keep1Candidates55 = chooseKiller1(5, 5, candidates55)

---

(word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1010110100, word: 1011100100, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1011010100, word: 1011100100, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1010110100, word: 1011100100, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1011010100, word: 1011100100, word: 1101110000)
0
this is the killer

---

(word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000)
0
this is the killer

============================================== 46 ==============================================

sage: keep1Candidates46 = chooseKiller1(4, 6, candidates46)

---

(word: 1101011010, word: 1101101010, word: 1110110100, word: 1110111000)
(1/24*u^4-5/12*u^3+35/24*u^2-25/12*u+1)*e[3, 3] + (1/120*u^5-1/8*u^4+17/24*u^3-15/8*u^2+137/60*u-1)*e[4, 2] + (1/720*u^6-1/48*u^5+17/144*u^4-5/16*u^3+137/360*u^2-1/6*u)*e[5, 1] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[6]


============================================== 56 ==============================================

sage: keep1Candidates56 = chooseKiller1(5, 6, candidates56)

---

(word: 11010101010, word: 11011010010, word: 11011101000, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11101010100, word: 11101011000, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11101010100, word: 11101101000, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11010110010, word: 11011100010, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11011010010, word: 11011100010, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11010110010, word: 11101011000, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11011010010, word: 11011101000, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11011001010, word: 11011100100, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]

---

(word: 11010101010, word: 11011001010, word: 11011100100, word: 11101110000)
(1/720*u^6-1/80*u^5+5/144*u^4-1/48*u^3-13/360*u^2+1/30*u)*e[3, 3] + (1/5040*u^7-1/360*u^6+1/72*u^5-1/36*u^4+7/720*u^3+11/360*u^2-1/42*u)*e[4, 2] + (1/40320*u^8-1/3360*u^7+1/960*u^6-11/1920*u^4+1/160*u^3+47/10080*u^2-1/168*u)*e[5, 1] + (1/362880*u^9-1/40320*u^8+1/60480*u^7+1/2880*u^6-11/17280*u^5-7/5760*u^4+59/22680*u^3+1/1120*u^2-1/504*u)*e[6]


============================================== 65 ==============================================

sage: keep1Candidates65 = chooseKiller1(6, 5, candidates65)
...

"""

# Return the polynomial expression, according to the shape of the last path and the maximal Hopf distance between the first two paths.
# This is the anklette statistic.
def refinedCountValidChainsPolynForbidden_r_q_first(m, n, forbiddenChains):
    return add([coeff * (binomial(r-2, length-1) + q^distances[0] * binomial(r-2, length-2)) * e(upStepsPartition(path)) for ((path, length, distances), coeff) in refinedCountValidChainsFirstLastForbidden(m, n, forbiddenChains).items()])

# Return the polynomial expression, according to the shape of the last path and the maximal Hopf distance between the last two paths.
# This is the collar statistic.
def refinedCountValidChainsPolynForbidden_r_q_last(m, n, forbiddenChains):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * e(upStepsPartition(path)) for ((path, length, distances), coeff) in refinedCountValidChainsFirstLastForbidden(m, n, forbiddenChains).items()])

# choose among the candidates using the statistics anklet and collar
def chooseKiller2(m, n, candidates):
    goal = e(Eval1(Phi[(m,n)], r-1+q, {r}))
    keepCandidates = []
    for candidate in candidates:
        print
        print "---"
        print
        print candidate
        forbiddenChains = strictSuperChains(TamariLattice(m,n), candidate, element_constructor=tuple)
        statsFirst = refinedCountValidChainsPolynForbidden_r_q_first(m,n, forbiddenChains)
        print statsFirst - goal
        statsLast = refinedCountValidChainsPolynForbidden_r_q_last(m,n, forbiddenChains)
        print statsLast - goal
        if statsFirst == goal and statsLast == goal:
            keepCandidates.append(candidate)
            print "this is the killer"
    return keepCandidates

"""

============================================== 55 ==============================================

sage: keep2Killers55 = chooseKiller2(5, 5, candidates55)

---

[word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000]
(-1/2*r^2*q^4+1/2*r^2*q^3+3/2*r*q^4-3/2*r*q^3-q^4+q^3)*e[3, 2] + (-1/6*r^3*q^4+1/6*r^3*q^3+r^2*q^4-r^2*q^3-11/6*r*q^4+11/6*r*q^3+q^4-q^3)*e[4, 1] + (-1/24*r^4*q^4+1/24*r^4*q^3+1/4*r^3*q^4-1/4*r^3*q^3-11/24*r^2*q^4+11/24*r^2*q^3+1/4*r*q^4-1/4*r*q^3)*e[5]
0

---

[word: 1010101010, word: 1010110100, word: 1011100100, word: 1101110000]
(-1/2*r^2*q^4+1/2*r^2*q^3+3/2*r*q^4-3/2*r*q^3-q^4+q^3)*e[3, 2] + (-1/6*r^3*q^4+1/6*r^3*q^3+r^2*q^4-r^2*q^3-11/6*r*q^4+11/6*r*q^3+q^4-q^3)*e[4, 1] + (-1/24*r^4*q^4+1/24*r^4*q^3+1/4*r^3*q^4-1/4*r^3*q^3-11/24*r^2*q^4+11/24*r^2*q^3+1/4*r*q^4-1/4*r*q^3)*e[5]
0

---

[word: 1010101010, word: 1011010100, word: 1011100100, word: 1101110000]
(-1/2*r^2*q^4+1/2*r^2*q^3+3/2*r*q^4-3/2*r*q^3-q^4+q^3)*e[3, 2] + (-1/6*r^3*q^4+1/6*r^3*q^3+r^2*q^4-r^2*q^3-11/6*r*q^4+11/6*r*q^3+q^4-q^3)*e[4, 1] + (-1/24*r^4*q^4+1/24*r^4*q^3+1/4*r^3*q^4-1/4*r^3*q^3-11/24*r^2*q^4+11/24*r^2*q^3+1/4*r*q^4-1/4*r*q^3)*e[5]
0

---

[word: 1010101010, word: 1011010100, word: 1101100100, word: 1110110000]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-7/12*r^3*q^3-3/2*r^2*q^4-r*q^5+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+71/24*r^2*q^3+16/3*r*q^4+2*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-77/12*r*q^3-6*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+5*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-7/240*r^5*q^2-1/6*r^4*q^3-7/12*r^3*q^4-1/2*r^2*q^5+11/2880*r^6+5/144*r^5*q+35/144*r^4*q^2+31/24*r^3*q^3+83/24*r^2*q^4+5/2*r*q^5-1/72*r^5-7/48*r^4*q-49/48*r^3*q^2-29/6*r^2*q^3-107/12*r*q^4-3*q^5+127/5760*r^4+29/90*r^3*q+203/90*r^2*q^2+87/10*r*q^3+8*q^4-1/288*r^3-7/20*r^2*q-49/20*r*q^2-6*q^3-29/1120*r^2+1/7*r*q+q^2+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-7/240*r^5*q^3-1/8*r^4*q^4-1/6*r^3*q^5+19/60480*r^7+11/2880*r^6*q+5/144*r^5*q^2+35/144*r^4*q^3+7/8*r^3*q^4+r^2*q^5-1/1440*r^6-1/72*r^5*q-7/48*r^4*q^2-49/48*r^3*q^3-23/8*r^2*q^4-11/6*r*q^5-11/17280*r^5+127/5760*r^4*q+29/90*r^3*q^2+203/90*r^2*q^3+247/60*r*q^4+q^5+13/2880*r^4-1/288*r^3*q-7/20*r^2*q^2-49/20*r*q^3-2*q^4-331/90720*r^3-29/1120*r^2*q+1/7*r*q^2+q^3-19/5040*r^2+1/56*r*q+1/252*r)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/360*r^6-1/48*r^5*q-1/12*r^4*q^2-1/4*r^3*q^3+1/72*r^5+11/144*r^4*q+7/24*r^3*q^2+11/24*r^2*q^3-1/36*r^4-1/16*r^3*q-5/12*r^2*q^2-1/4*r*q^3+7/720*r^3-7/90*r^2*q+1/5*r*q^2+11/360*r^2+1/12*r*q-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-1/48*r^5*q^2-1/8*r^4*q^3-1/4*r^3*q^4+11/2880*r^6+19/720*r^5*q+17/144*r^4*q^2+13/24*r^3*q^3+11/24*r^2*q^4-1/72*r^5-1/16*r^4*q-5/16*r^3*q^2-7/8*r^2*q^3-1/4*r*q^4+127/5760*r^4+11/360*r^3*q+137/360*r^2*q^2+9/20*r*q^3-1/288*r^3+1/15*r^2*q-1/6*r*q^2-29/1120*r^2-2/35*r*q+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5+19/60480*r^7+11/2880*r^6*q+19/720*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5-1/1440*r^6-1/72*r^5*q-1/16*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5-11/17280*r^5+127/5760*r^4*q+11/360*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4+13/2880*r^4-1/288*r^3*q+1/15*r^2*q^2-1/6*r*q^3-331/90720*r^3-29/1120*r^2*q-2/35*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]

---

[word: 1010101010, word: 1101010100, word: 1101101000, word: 1110110000]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-5/12*r^3*q^3-r^2*q^4+1/72*r^5+17/144*r^4*q+13/24*r^3*q^2+35/24*r^2*q^3+11/6*r*q^4-1/36*r^4-5/16*r^3*q-7/8*r^2*q^2-25/12*r*q^3-q^4+7/720*r^3+137/360*r^2*q+9/20*r*q^2+q^3+11/360*r^2-1/6*r*q-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-7/240*r^5*q^2-1/8*r^4*q^3-5/12*r^3*q^4+11/2880*r^6+5/144*r^5*q+29/144*r^4*q^2+17/24*r^3*q^3+35/24*r^2*q^4-1/72*r^5-7/48*r^4*q-29/48*r^3*q^2-15/8*r^2*q^3-25/12*r*q^4+127/5760*r^4+29/90*r^3*q+287/360*r^2*q^2+137/60*r*q^3+q^4-1/288*r^3-7/20*r^2*q-11/30*r*q^2-q^3-29/1120*r^2+1/7*r*q+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-1/48*r^5*q^3-1/12*r^4*q^4+19/60480*r^7+11/2880*r^6*q+19/720*r^5*q^2+17/144*r^4*q^3+7/24*r^3*q^4-1/1440*r^6-1/72*r^5*q-1/16*r^4*q^2-5/16*r^3*q^3-5/12*r^2*q^4-11/17280*r^5+127/5760*r^4*q+11/360*r^3*q^2+137/360*r^2*q^3+1/5*r*q^4+13/2880*r^4-1/288*r^3*q+1/15*r^2*q^2-1/6*r*q^3-331/90720*r^3-29/1120*r^2*q-2/35*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-5/12*r^3*q^3-1/2*r^2*q^4+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+47/24*r^2*q^3+3/2*r*q^4-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-43/12*r*q^3-q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+2*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-7/240*r^5*q^2-1/6*r^4*q^3-5/12*r^3*q^4-1/2*r^2*q^5+11/2880*r^6+5/144*r^5*q+35/144*r^4*q^2+23/24*r^3*q^3+47/24*r^2*q^4+3/2*r*q^5-1/72*r^5-7/48*r^4*q-41/48*r^3*q^2-7/3*r^2*q^3-43/12*r*q^4-q^5+127/5760*r^4+29/90*r^3*q+113/90*r^2*q^2+38/15*r*q^3+2*q^4-1/288*r^3-7/20*r^2*q-37/60*r*q^2-q^3-29/1120*r^2+1/7*r*q+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-7/240*r^5*q^3-1/6*r^4*q^4-5/12*r^3*q^5-1/2*r^2*q^6+19/60480*r^7+11/2880*r^6*q+5/144*r^5*q^2+35/144*r^4*q^3+23/24*r^3*q^4+47/24*r^2*q^5+3/2*r*q^6-1/1440*r^6-1/72*r^5*q-7/48*r^4*q^2-41/48*r^3*q^3-7/3*r^2*q^4-43/12*r*q^5-q^6-11/17280*r^5+127/5760*r^4*q+29/90*r^3*q^2+113/90*r^2*q^3+38/15*r*q^4+2*q^5+13/2880*r^4-1/288*r^3*q-7/20*r^2*q^2-37/60*r*q^3-q^4-331/90720*r^3-29/1120*r^2*q+1/7*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]

---

[word: 1010101010, word: 1010110100, word: 1011011000, word: 1110011000]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-7/12*r^3*q^3-3/2*r^2*q^4-r*q^5-q^6+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+71/24*r^2*q^3+13/3*r*q^4+3*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-65/12*r*q^3-4*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+3*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-7/240*r^5*q^2-1/6*r^4*q^3-7/12*r^3*q^4-1/2*r^2*q^5-r*q^6+11/2880*r^6+5/144*r^5*q+35/144*r^4*q^2+31/24*r^3*q^3+71/24*r^2*q^4+7/2*r*q^5+2*q^6-1/72*r^5-7/48*r^4*q-49/48*r^3*q^2-13/3*r^2*q^3-77/12*r*q^4-5*q^5+127/5760*r^4+29/90*r^3*q+203/90*r^2*q^2+31/5*r*q^3+5*q^4-1/288*r^3-7/20*r^2*q-49/20*r*q^2-3*q^3-29/1120*r^2+1/7*r*q+q^2+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-7/240*r^5*q^3-1/8*r^4*q^4-1/6*r^3*q^5-1/2*r^2*q^6+19/60480*r^7+11/2880*r^6*q+5/144*r^5*q^2+35/144*r^4*q^3+17/24*r^3*q^4+3/2*r^2*q^5+3/2*r*q^6-1/1440*r^6-1/72*r^5*q-7/48*r^4*q^2-41/48*r^3*q^3-15/8*r^2*q^4-10/3*r*q^5-q^6-11/17280*r^5+127/5760*r^4*q+29/90*r^3*q^2+113/90*r^2*q^3+137/60*r*q^4+2*q^5+13/2880*r^4-1/288*r^3*q-7/20*r^2*q^2-37/60*r*q^3-q^4-331/90720*r^3-29/1120*r^2*q+1/7*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-5/12*r^3*q^3-1/2*r^2*q^4-r*q^5-q^6+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+47/24*r^2*q^3+5/2*r*q^4+3*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-43/12*r*q^3-3*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+2*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-1/48*r^5*q^2-1/8*r^4*q^3-1/4*r^3*q^4+11/2880*r^6+19/720*r^5*q+17/144*r^4*q^2+13/24*r^3*q^3+11/24*r^2*q^4-1/72*r^5-1/16*r^4*q-5/16*r^3*q^2-7/8*r^2*q^3-1/4*r*q^4+127/5760*r^4+11/360*r^3*q+137/360*r^2*q^2+9/20*r*q^3-1/288*r^3+1/15*r^2*q-1/6*r*q^2-29/1120*r^2-2/35*r*q+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5+19/60480*r^7+11/2880*r^6*q+19/720*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5-1/1440*r^6-1/72*r^5*q-1/16*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5-11/17280*r^5+127/5760*r^4*q+11/360*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4+13/2880*r^4-1/288*r^3*q+1/15*r^2*q^2-1/6*r*q^3-331/90720*r^3-29/1120*r^2*q-2/35*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]

---

[word: 1010101010, word: 1010110100, word: 1101011000, word: 1110011000]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3+1/6*r^3*q^4-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-7/12*r^3*q^3-3/2*r^2*q^4-r*q^5-q^6+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+71/24*r^2*q^3+16/3*r*q^4+3*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-77/12*r*q^3-6*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+5*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-7/240*r^5*q^2-1/6*r^4*q^3-7/12*r^3*q^4-1/2*r^2*q^5-r*q^6+11/2880*r^6+5/144*r^5*q+35/144*r^4*q^2+31/24*r^3*q^3+83/24*r^2*q^4+7/2*r*q^5+2*q^6-1/72*r^5-7/48*r^4*q-49/48*r^3*q^2-29/6*r^2*q^3-107/12*r*q^4-5*q^5+127/5760*r^4+29/90*r^3*q+203/90*r^2*q^2+87/10*r*q^3+8*q^4-1/288*r^3-7/20*r^2*q-49/20*r*q^2-6*q^3-29/1120*r^2+1/7*r*q+q^2+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-7/240*r^5*q^3-1/8*r^4*q^4-1/6*r^3*q^5-1/2*r^2*q^6+19/60480*r^7+11/2880*r^6*q+5/144*r^5*q^2+35/144*r^4*q^3+7/8*r^3*q^4+3/2*r^2*q^5+3/2*r*q^6-1/1440*r^6-1/72*r^5*q-7/48*r^4*q^2-49/48*r^3*q^3-23/8*r^2*q^4-10/3*r*q^5-q^6-11/17280*r^5+127/5760*r^4*q+29/90*r^3*q^2+203/90*r^2*q^3+247/60*r*q^4+2*q^5+13/2880*r^4-1/288*r^3*q-7/20*r^2*q^2-49/20*r*q^3-2*q^4-331/90720*r^3-29/1120*r^2*q+1/7*r*q^2+q^3-19/5040*r^2+1/56*r*q+1/252*r)*e[5]
(1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/360*r^6-1/48*r^5*q-1/8*r^4*q^2-5/12*r^3*q^3-1/2*r^2*q^4-r*q^5-q^6+1/72*r^5+17/144*r^4*q+17/24*r^3*q^2+47/24*r^2*q^3+5/2*r*q^4+3*q^5-1/36*r^4-5/16*r^3*q-15/8*r^2*q^2-43/12*r*q^3-3*q^4+7/720*r^3+137/360*r^2*q+137/60*r*q^2+2*q^3+11/360*r^2-1/6*r*q-q^2-1/42*r)*e[3, 2] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/2016*r^7-1/240*r^6*q-1/48*r^5*q^2-1/8*r^4*q^3-1/4*r^3*q^4+11/2880*r^6+19/720*r^5*q+17/144*r^4*q^2+13/24*r^3*q^3+11/24*r^2*q^4-1/72*r^5-1/16*r^4*q-5/16*r^3*q^2-7/8*r^2*q^3-1/4*r*q^4+127/5760*r^4+11/360*r^3*q+137/360*r^2*q^2+9/20*r*q^3-1/288*r^3+1/15*r^2*q-1/6*r*q^2-29/1120*r^2-2/35*r*q+1/56*r)*e[4, 1] + (1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5+19/60480*r^7+11/2880*r^6*q+19/720*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5-1/1440*r^6-1/72*r^5*q-1/16*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5-11/17280*r^5+127/5760*r^4*q+11/360*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4+13/2880*r^4-1/288*r^3*q+1/15*r^2*q^2-1/6*r*q^3-331/90720*r^3-29/1120*r^2*q-2/35*r*q^2-19/5040*r^2+1/56*r*q+1/252*r)*e[5]

============================================== 46 ==============================================

sage: keep2Killers46 = chooseKiller2(4, 6, keep1Candidates46)
NOT UPDATED
---

[word: 1101011010, word: 1101011100, word: 1110101100, word: 1110111000]
0
this is the killer

---

[word: 1101011010, word: 1101101100, word: 1110110100, word: 1110111000]
(1/2*u^2*v^2-1/2*u^2*v-5/2*u*v^2+5/2*u*v+3*v^2-3*v)*e[3, 3] + (1/2*u^2*v^3-u^2*v^2-5/2*u*v^3+1/2*u^2*v+5*u*v^2+3*v^3-5/2*u*v-6*v^2+3*v)*e[4, 2] + (1/2*u^2*v^4-u^2*v^3-5/2*u*v^4+1/2*u^2*v^2+5*u*v^3+3*v^4-5/2*u*v^2-6*v^3+3*v^2)*e[5, 1] + (1/2*u^2*v^5-u^2*v^4-5/2*u*v^5+1/2*u^2*v^3+5*u*v^4+3*v^5-5/2*u*v^3-6*v^4+3*v^3)*e[6]

---

[word: 1101011010, word: 1110101100, word: 1110110100, word: 1110111000]
(-u*v^2+u*v+2*v^2-2*v)*e[3, 3] + (-1/2*u^2*v^2-u*v^3+1/2*u^2*v+7/2*u*v^2+2*v^3-5/2*u*v-5*v^2+3*v)*e[4, 2] + (-1/2*u^2*v^3-u*v^4+1/2*u^2*v^2+7/2*u*v^3+2*v^4-5/2*u*v^2-5*v^3+3*v^2)*e[5, 1] + (-1/2*u^2*v^4-u*v^5+1/2*u^2*v^3+7/2*u*v^4+2*v^5-5/2*u*v^3-5*v^4+3*v^3)*e[6]

---

[word: 1101011010, word: 1101101010, word: 1101110010, word: 1110111000]
0
this is the killer

---

[word: 1101011010, word: 1101101010, word: 1101110100, word: 1110111000]
0
this is the killer

---

[word: 1101011010, word: 1101110010, word: 1101111000, word: 1110111000]
(-u*v^2+u*v+2*v^2-2*v)*e[3, 3] + (-1/2*u^2*v^2-u*v^3+1/2*u^2*v+7/2*u*v^2+2*v^3-5/2*u*v-5*v^2+3*v)*e[4, 2] + (-1/2*u^2*v^3-u*v^4+1/2*u^2*v^2+7/2*u*v^3+2*v^4-5/2*u*v^2-5*v^3+3*v^2)*e[5, 1] + (-1/2*u^2*v^4-u*v^5+1/2*u^2*v^3+7/2*u*v^4+2*v^5-5/2*u*v^3-5*v^4+3*v^3)*e[6]


============================================== 56 ==============================================

sage: keep2Killers56 = chooseKiller2(5, 6, keep1Candidates56)
NOT UPDATED
---

[word: 11010101010, word: 11010110100, word: 11101011000, word: 11101110000]
0
this is the killer

---

[word: 11010101010, word: 11010110100, word: 11011100100, word: 11101110000]
0
this is the killer

---

[word: 11010101010, word: 11011010100, word: 11011100100, word: 11101110000]
0
this is the killer


============================================== 65 ==============================================

sage: keep2Killers65 = chooseKiller2(6, 5, keep1Candidates65)
NOT UPDATED
---

[word: 10101010100, word: 10101101000, word: 11010110000, word: 11011100000]
0
this is the killer

---

[word: 10101010100, word: 10101101000, word: 10111001000, word: 11011100000]
0
this is the killer

---

[word: 10101010100, word: 10110101000, word: 10111001000, word: 11011100000]
0
this is the killer

"""

# returns the polynomial expression
def refinedCountValidChainsPolynForbidden_r_q_LLT_last(m, n, forbiddenChains):
    return add([coeff * (binomial(r-2, length-1) + q^distances[-1] * binomial(r-2, length-2)) * LLT[partitionPath(path)[:-1]] for ((path, length, distances), coeff) in refinedCountValidChainsFirstLastForbidden(m, n, forbiddenChains).items()])

# choose among the candidates using the statistic
def chooseKiller3(m, n, candidates):
    goal = Eval1(Phi[(n,n)], q+t+r-2, {r})
    keepCandidates = []
    for candidate in candidates:
        print
        print "---"
        print
        print candidate
        forbiddenChains = strictSuperChains(TamariLattice(m,n), candidate, element_constructor=tuple)
        statsLast = refinedCountValidChainsPolynForbidden_r_q_LLT_last(m,n, forbiddenChains)
        print statsLast - goal
        if statsLast == goal:
            keepCandidates.append(candidate)
            print "this is the killer"
    return keepCandidates

"""
sage: keep3Killers55 = chooseKiller3(5, 5, candidates55)

---

[word: 1010101010, word: 1010110100, word: 1101011000, word: 1101110000]
0
this is the killer

---

[word: 1010101010, word: 1010110100, word: 1011100100, word: 1101110000]
0
this is the killer

---

[word: 1010101010, word: 1011010100, word: 1011100100, word: 1101110000]
0
this is the killer

---

[word: 1010101010, word: 1011010100, word: 1101100100, word: 1110110000]
(1/362880*r^9+1/40320*r^8*q+1/5040*r^7*q^2+1/720*r^6*q^3+1/120*r^5*q^4+1/24*r^4*q^5+1/40320*r^8*t+1/5040*r^7*q*t+1/720*r^6*q^2*t+1/120*r^5*q^3*t+1/24*r^4*q^4*t-1/20160*r^8-1/2016*r^7*q-1/240*r^6*q^2-1/48*r^5*q^3-1/8*r^4*q^4-1/4*r^3*q^5-1/3360*r^7*t-1/360*r^6*q*t-1/80*r^5*q^2*t-1/12*r^4*q^3*t-1/4*r^3*q^4*t+19/60480*r^7+11/2880*r^6*q+19/720*r^5*q^2+17/144*r^4*q^3+13/24*r^3*q^4+11/24*r^2*q^5+1/960*r^6*t+1/180*r^5*q*t+5/144*r^4*q^2*t+7/24*r^3*q^3*t+11/24*r^2*q^4*t-1/1440*r^6-1/72*r^5*q-1/16*r^4*q^2-5/16*r^3*q^3-7/8*r^2*q^4-1/4*r*q^5+1/72*r^4*q*t-1/48*r^3*q^2*t-5/12*r^2*q^3*t-1/4*r*q^4*t-11/17280*r^5+127/5760*r^4*q+11/360*r^3*q^2+137/360*r^2*q^3+9/20*r*q^4-11/1920*r^4*t-23/720*r^3*q*t-13/360*r^2*q^2*t+1/5*r*q^3*t+13/2880*r^4-1/288*r^3*q+1/15*r^2*q^2-1/6*r*q^3+1/160*r^3*t-1/90*r^2*q*t+1/30*r*q^2*t-331/90720*r^3-29/1120*r^2*q-2/35*r*q^2+47/10080*r^2*t+11/420*r*q*t-19/5040*r^2+1/56*r*q-1/168*r*t+1/252*r)*s[1, 1, 1, 1, 1] + (1/40320*r^8+1/5040*r^7*q+1/720*r^6*q^2+1/120*r^5*q^3+1/24*r^4*q^4-1/3360*r^7-1/360*r^6*q-1/80*r^5*q^2-1/12*r^4*q^3-1/4*r^3*q^4+1/960*r^6+1/180*r^5*q+5/144*r^4*q^2+7/24*r^3*q^3+11/24*r^2*q^4+1/72*r^4*q-1/48*r^3*q^2-5/12*r^2*q^3-1/4*r*q^4-11/1920*r^4-23/720*r^3*q-13/360*r^2*q^2+1/5*r*q^3+1/160*r^3-1/90*r^2*q+1/30*r*q^2+47/10080*r^2+11/420*r*q-1/168*r)*s[2, 1, 1, 1] + (1/5040*r^7+1/720*r^6*q+1/120*r^5*q^2+1/24*r^4*q^3-1/360*r^6-1/48*r^5*q-1/12*r^4*q^2-1/4*r^3*q^3+1/72*r^5+11/144*r^4*q+7/24*r^3*q^2+11/24*r^2*q^3-1/36*r^4-1/16*r^3*q-5/12*r^2*q^2-1/4*r*q^3+7/720*r^3-7/90*r^2*q+1/5*r*q^2+11/360*r^2+1/12*r*q-1/42*r)*s[2, 2, 1]

---

ALL OTHER ARE BAD
"""

""" ======================= LOOKING FOR A KILLER FAMILY ========================= """

"""
    Now we try to understand the possible killer families:
    - compute the chains that are killed by a given family of killers
    - find candidates for the killers
"""

def candidateKillerFamilies(m, n, mu):
    killerFamilies = []
    bp = bottomPath(m,n) # the bottom path
    # first, we look for elements that kill a given excess
    killerDictionary = defaultdict(lambda: tuple([]))
    for tp in upStepDictionnary(m,n)[mu]: # the possible top paths
        belowtp = TamariLattice(m,n).subposet(TamariLattice(m,n).principal_lower_set(tp))
        for j in range(4,7):
            for candidateKiller in [tuple(c) for c in belowtp.chains() if len(c) == j and c[0] == bp and c[-1] == tp]:
                killedChains = tuple(strictSuperChains(belowtp, candidateKiller, element_constructor=tuple))
                contribution = sum(x^len(c) for c in killedChains if isValidChain(c))
                for i in range(1,8-j):
                    goal = add(x^len(c) for c in Subsets(range(i))) * x^j
                    if contribution == goal:
                        killerDictionary[(i,j)] += ((candidateKiller, killedChains),)
    # a short print to see how difficult it will be...
    print 'here is the size of the dictionary...'
    for i in range(1,4):
        for j in range(4,5+i):
            print i, j, len(killerDictionary[(i,j)])
    print '...'
    # second, we try to combine killers together
    for p in [[(2,4),(2,5)]]:
        print p
        for candidateFamily in cartesian_product([killerDictionary[ij] for ij in p]):
            killers = [candidateKiller for (candidateKiller, killedChains) in candidateFamily]
            killedChains = reduce(union, [killedChains for (candidateKiller, killedChains) in candidateFamily])
            contribution = sum(x^len(c) for c in killedChains if isValidChain(c))
            if contribution == add(x^len(c) for c in Subsets(range(3))) * x^4:
                # print candidateFamily
                killerFamilies.append((killers, killedChains))
    return killerFamilies

# choose among the candidate families using the statistics anklet and collar
def chooseKillerFamily2(m, n, candidateKillerFamilies):
    goal = e(Eval1(Phi[(m,n)], r-1+q, {r}))
    keepCandidates = []
    cmp = 0
    for (candidateFamily, killedChains) in candidateKillerFamilies:
        print
        print "---"
        print cmp
        print candidateFamily
        statsFirst = refinedCountValidChainsPolynForbidden_r_q_first(m,n, killedChains)
        print statsFirst - goal
        statsLast = refinedCountValidChainsPolynForbidden_r_q_last(m,n, killedChains)
        print statsLast - goal
        if statsFirst == goal and statsLast == goal:
            keepCandidates.append(candidateFamily)
            print "this is the killer"
        cmp += 1
    return keepCandidates

