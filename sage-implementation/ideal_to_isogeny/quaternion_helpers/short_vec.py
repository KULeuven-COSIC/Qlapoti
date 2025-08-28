from sage.all import *
import fpylll
from fpylll import IntegerMatrix, CVP
from fpylll.fplll.gso import MatGSO


#####################################################
#                                                   #
#    Enumerate close vectors, from LearningToSQI    #
# https://github.com/LearningToSQI/SQISign-SageMath #
#                                                   #
#####################################################

"""MIT License

Copyright (c) 2023 Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer  and Giacomo Pope

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software."""


def solve_closest_vector_problem(lattice_basis, target, reduced=False):
    """
    Use the fpylll library to solve the CVP problem for a given
    lattice basis and target vector
    """
    if reduced:
        L = IntegerMatrix.from_matrix(lattice_basis)
    else:
        L = IntegerMatrix.from_matrix(lattice_basis.LLL())
    v = CVP.closest_vector(L, target)
    # fpylll returns a type `tuple` object
    return vector(v)


def generate_short_vectors_fpyll(L, bound, count=2000):
    """
    Helper function for GenerateShortVectors and 
    generate_small_norm_quat which builds an iterator
    for short norm vectors of an LLL reduced lattice
    basis.
    """
    # # Move from Sage world to Fypll world
    A = IntegerMatrix.from_matrix(L)

    # Gram-Schmidt Othogonalization
    G = MatGSO(A)
    _ = G.update_gso()

    # Enumeration class on G with `count`` solutions
    # BEST_N_SOLUTIONS:
    # Starting with the nr_solutions-th solution, every time a new solution is found
    # the enumeration bound is updated to the length of the longest solution. If more
    # than nr_solutions were found, the longest is dropped.
    E = fpylll.Enumeration(
        G, nr_solutions=count, strategy=fpylll.EvaluatorStrategy.BEST_N_SOLUTIONS
    )

    # We need the row count when we call enumerate
    r = L.nrows()

    # If enumerate finds no solutions it raises an error, so we
    # wrap it in a try block
    try:
        # The arguments of enumerate are:
        # E.enumerate(first_row, last_row, max_dist, max_dist_expo)
        short_vectors = E.enumerate(0, r, bound, 0)
    except Exception as e:
        short_vectors = []
        
    return short_vectors

def generate_short_vectors(lattice_basis, bound, reduced=False, count=2000):
    """
    Generate a generator of short vectors with norm <= `bound`
    returns at most `count` vectors.
    
    Most of the heavy lifting of this function is done by 
    generate_short_vectors_fpyll
    """
    if reduced:
        L = lattice_basis
    else:
        L = lattice_basis.LLL()
    
    short_vectors = generate_short_vectors_fpyll(L, bound, count=count)
    for _, xis in short_vectors:
        # Returns values x1,x2,...xr such that
        # x0*row[0] + ... + xr*row[r] = short vector
        v3 = vector([ZZ(xi) for xi in xis])
        v = v3 * L
        yield v


def generate_close_vectors(lattice_basis, target, p, L, reduced=False, count=1000):
    """
    Generate a generator of vectors which are close, without
    bound determined by N to the `target`. The first
    element of the list is the solution of the CVP.
    """
    # Compute the closest element
    closest = solve_closest_vector_problem(lattice_basis, target, reduced=reduced)
    yield closest

    # Now use short vectors below a bound to find
    # close enough vectors

    # Set the distance
    diff = target - closest
    distance = diff.dot_product(diff)

    # Compute the bound from L
    #I don't understand why the bound should get bigger if the shortest vector is very close....
    b0 = L // p
    #bound = floor((b0 + distance) + (2 * (b0 * distance).sqrt()))
    bound = floor(b0)

    short_vectors = generate_short_vectors(lattice_basis, bound, reduced=reduced, count=count)

    counter = 0
    for v in short_vectors:
        if counter > 1000:
            print("????")
            break
        counter += 1
        yield closest + v