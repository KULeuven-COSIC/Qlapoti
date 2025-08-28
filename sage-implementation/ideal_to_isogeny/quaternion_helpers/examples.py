from helpers import *

def test_successive_minima():
    p = 10007
    B = QuaternionAlgebra(-1, -p)
    O0 = B.maximal_order()

    #Generate random orders, compute their successive minima
    for _ in range(10):
        I = heuristic_random_ideal(O0, next_prime(randint(p, p**2)))
        O = I.right_order()

        mins = successive_minima(O)
        print("\n~~~~~~~~~~~~~~~~")
        print(O)
        print(f"successive minima: {mins}")
        print(f"successive minima log_p: {[round(log(lam, p), 5) for lam in mins]}")

def test_basis_thing():
    p = 10007
    B = QuaternionAlgebra(-1, -p)
    O0 = B.maximal_order()
    ell = next_prime(randint(p, p**2))
    I = heuristic_random_ideal(O0, ell)
    alpha = ideal_generator(I)
    gamma = O0.random_element()*alpha + O0.random_element()*ell
    beta_1, beta_2 = decompose_in_ideal(gamma, alpha, ell, O0)
    assert beta_1*alpha + beta_2*ell == gamma
    print(f"success!")

def test_pushforward_endo_thing():
    p = 10007
    B = QuaternionAlgebra(-1, -p)
    O0 = B.maximal_order()
    ell = next_prime(randint(p, p**2))

    print("Testing two random ideals")
    I = heuristic_random_ideal(O0, ell)
    J = heuristic_random_ideal(O0, ell)

    theta = find_correct_pushforward(I, J)
    J_prime = pushforward(I, (O0*theta))

    assert is_isomorphic(J_prime.left_order(), J.left_order())
    assert is_isomorphic(J_prime.right_order(), J.right_order())

    print("     > Success!")
    print("Testing where one is an endomorphism")
    rho = represent_integer(O0, ell)
    J = O0*rho

    theta = find_correct_pushforward_endo(I, J)
    J_prime = pushforward(I, (O0*theta))

    assert is_isomorphic(J_prime.left_order(), J.left_order())
    assert is_isomorphic(J_prime.right_order(), J.right_order())
    
    print("     > Success!")

def test_pushforward_thing():
    p = 10007
    B = QuaternionAlgebra(-1, -p)
    O0 = B.maximal_order()
    ell = next_prime(randint(p, p**2))

    print("Testing two random ideals")
    I = heuristic_random_ideal(O0, ell).conjugate()
    J = O0*represent_integer(O0, ell)

    K = find_correct_pushforward(I, J)
    J_prime = pushforward(I, K)

    assert is_isomorphic(J_prime.left_order(), J.left_order())
    assert is_isomorphic(J_prime.right_order(), J.right_order())

    print("     > Success!")

    

if __name__ == "__main__":
    test_pushforward_thing()
    