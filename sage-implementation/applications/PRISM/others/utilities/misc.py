def randomise_basis():
    # generating a torsion basis
    P0_, Q0_ = torsion_basis_2e(E0, e_prime)
    # generate a random basis
    a = randint(0, 2**e_prime)
    b = randint(0, 2**e_prime)
    c = randint(0, 2**e_prime)
    d = randint(0, 2**e_prime)

    while (a * d - b * c) % 2 == 0:
        a = randint(0, 2**e_prime)
        b = randint(0, 2**e_prime)
        c = randint(0, 2**e_prime)
        d = randint(0, 2**e_prime)

    P0 = a * P0_ + b * Q0_
    Q0 = c * P0_ + d * Q0_

    return P0, Q0
