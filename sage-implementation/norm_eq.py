from sage.all import *
from ideal_to_isogeny.quaternion_helpers.helpers import *
from applications.PRISM.utilities.sum_of_squares import sum_of_squares
from time import time

def smallish_gen(I, N, I_basis=None):
    if not I_basis:
        print("SHOULDNT GET CALLED")
        I_basis = reduced_basis(I)
    #N = I.norm() This is slow, pass instead

    while True:
        alpha = sum(randint(1,100000)*gen for gen in I_basis)

        a_alpha = alpha.coefficient_tuple()[0]
        if gcd(2*a_alpha, N) == 1 and gcd(alpha.reduced_norm(), N*N) == N:
            break

    assert gcd(alpha.reduced_norm(), N**2) == N
    assert alpha in I
    return alpha

def _succ_min(L, p):
    fourth_root_p = round(p**(1/4), 10)
    lam1 = round(L.row(0).norm()/fourth_root_p, 10)
    lam2 = round(L.row(1).norm()/fourth_root_p, 10)
    return lam1, lam2

def predict_endtype(tlist):

    #transformlist = {1: [[2, 0, 1, 2], [2, 2, 3, 0], [2, 2, 1, 0], [2, 0, 3, 2]], 2: [[0, 2, 2, 1], [0, 2, 2, 3], [2, 2, 0, 1], [2, 2, 0, 3]]}
    # t1, t2, t3, t4 = tlist
    # Again, disagreeing on what j and k is
    t1, t2, t4, t3 = tlist

    if t1 == 2:
        if t2 == 2:
            if t3 == 1 and t4 == 0:
                return 1
            elif t3 == 3 and t4 == 2:
                return 1
            elif t3 == 0:
                if t4 == 1:
                    return 2
                elif t4 == 3:
                    return 2
                else:
                    return 0
            else:
                return 0
        elif t2 == 0:
            if t3 == 1 and t4 == 2:
                return 1
            elif t3 == 3 and t4 == 2:
                return 1
            else:
                return 0
        else:
            return 0
    elif t1 == 0:
        if t2 == 2 and t3 == 2:
            if t4 == 1:
                return 2
            elif t4 == 3:
                return 2
            else:
                return 0
        else:
            return 0
    else:
        return 0

def qlapoti(J, e, allow_diags = True, odd_norm_output = False, stats = False):
    I, betaij = reduced_ideal(J, return_elt=True) #The smallest ideal in the class
    p = I.quaternion_algebra().ramified_primes()[0]
    num_tries = 0
    num_cornacchia = 0

    # The algorithm...
    N = ZZ(I.norm())

    solved = False
    I_basis = reduced_basis(I)

    keep_alpha = False
    lam = 1


    while not solved:
        num_tries += 1
        if num_tries % 20000 == 0:
            lam = 1
            keep_alpha = False
        #print("\n\nTrying new gen....")
        if keep_alpha:
            alpha += alpha_0
            lam += 1
            a_alpha = list(alpha.coefficient_tuple())[0]
            #multiply alpha by an invertible lambda, and check that 2*alpha is still invertible
            while gcd(lam, N) > 1:
                alpha += alpha_0
                lam += 1
                a_alpha = list(alpha.coefficient_tuple())[0]
        else:
            alpha_0 = smallish_gen(I, N, I_basis=I_basis)
            alpha = alpha_0
            #print(f"alpha = {alpha}")
            alpha_0_norm = alpha.reduced_norm()
            #print(f"norm = {factor(ZZ(alpha_0_norm), limit=100)}")

        a_alpha = list(alpha.coefficient_tuple())[0]
        b_alpha = list(alpha.coefficient_tuple())[1]

        M = ZZ(2**e - 2*lam**2*alpha_0_norm/N)
        if M < 0:
            print("M already too small! NB! Shouldnt happen")
            continue
        #N(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha1*a1 + 2*a_alpha2*a2 + 2*balpha1*b1 + 2*balpha2*b2 = M
        #b1 + balpha1^-1*balpha2*b2) = M*(2*balpha1)^-1 (mod N)
        #A + (2*b_alpha/2*a_alpha)B = M/(2*a_alpha) (mod N)

        if keep_alpha:
            a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha), N)
            T = (M*a_alpha_inv2) % N
            v_target = vector(ZZ, [-T, 0])
        else:
            a_alpha_inv2 = inverse_mod(ZZ(2*a_alpha), N)
            x = ZZ((2*b_alpha*a_alpha_inv2) % N)
            L = Matrix(ZZ, [[N-x, 1], [N, 0]])
            L = L.LLL() # Make custom gauss reduction if necessary for speedup...

            T = ZZ((M*a_alpha_inv2) % N)
            v_target = vector(ZZ, [-T, 0])
            L_inv = Matrix(QQ, L).inverse()
            if _succ_min(L, p)[1] < 1 and gcd(2*a_alpha, 2*b_alpha) == 1: #Just make sure our lattice isnt completely unreasonable
                if not stats:
                    keep_alpha = True

        AB_vec = vector([round(c) for c in v_target*L_inv])
        v_close = AB_vec*L
    
        A, B = v_close - v_target

        assert (2*(a_alpha*A + b_alpha*B)) % N == M % N

        M2 = M - 2*a_alpha*A - 2*b_alpha*B
        M2 = ZZ(M2/N)

        # Complete the square
        # (a2 = A - a1, b2 = B - b1)
        #2*a1^2 + A^2 - 2Aa1 + 2*b1^2 + B^2 - 2Bb1 = M2
        #M3 = M2 - A**2 - B**2

        #(2*a1)^2 + (2*b1)^2 - 4*Aa1 - 4*Bb1 = (2*a1 - A)^2 + (2*b1 - B)^2 = 4*M3 + A^2 + B^2
        # 4*M3 = (2*a1)^2 + (2*b1)^2 - 4*A*a1 - 4*B*b1 = (2*a1 - A)^2 + (2*b1 - B)^2 - (A^2 + B^2)

        M4 = 2*M2 - A**2 - B**2
        if not allow_diags and M4 < 0:
            continue

        #Unsolvable cases
        if M4 % 8 == 0:
            continue
        if A % 2 == B % 2 == 0:
            if M4 % 4 != 0:
                continue
        elif A % 2 == B % 2 == 1:
            if M4 % 4 != 2:
                continue
        else:
            if M4 % 4 != 1:
                continue

        #### Can we predict already here?
        
        num_cornacchia += 1

        ab = sum_of_squares(M4)

        if not ab:
            continue

        solved = True
        ad1, bd1 = ab

        if ad1 % 2 != A % 2:
            temp = ad1
            ad1 = bd1
            bd1 = temp

        assert ad1 % 2 == A % 2 and bd1 % 2 == B % 2
        solved = True
        a1 = ZZ((ad1 + A)/2)
        b1 = ZZ((bd1 + B)/2)
        a2 = A - a1
        b2 = B - b1

        i, j, k = I.quaternion_algebra().gens()
        gamma1 = a1 + i*b1
        gamma2 = a2 + i*b2

        mu1 = N*gamma1 + alpha
        mu2 = N*gamma2 + alpha

        θ = mu2 * mu1.conjugate() / N

        if allow_diags:
            d1 = mu1.reduced_norm()/N
            d2 = mu2.reduced_norm()/N
            if not odd_norm_output or d1 % 2 == d2 % 2 == 1:
                break
            else:
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
        
        if θ.denominator() == 1:
            # If this should be possible, it must come from a "transformation"
            theta_profile = [c % 4 for c in θ]
            res = predict_endtype(theta_profile)
            if res == 0:
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
            elif res == 1:
                temp = (mu1 + mu2)/2
                mu2 = (mu1 - mu2)/2
                mu1 = temp
            else:
                temp = (mu1 + i*mu2)/2
                mu2 = (mu1 - i*mu2)/2
                mu1 = temp
        
            e = e-1 
            θ = mu2 * mu1.conjugate() / N
        else:
            # Check odd norm
            theta_profile = [c % 4 for c in 2*θ]
            tx, ty, tz, tw = 2*θ 
            if not ((tx % 2 != tz % 2) and ([c % 4 for c in 2*θ].count(2) == 1)):
                solved = False
                keep_alpha = False #Some congruence condition it seems like, so can get stuck
                lam = 1
                continue
                
        d1 = mu1.reduced_norm()/N
        d2 = mu2.reduced_norm()/N

    assert allow_diags or θ.denominator() == 2 
    
    assert solved
    
    assert N*(a1**2 + a2**2 + b1**2 + b2**2) + 2*a_alpha*(a1 + a2) + 2*b_alpha*(b1 + b2) == M

    assert mu1 in I
    assert mu2 in I
    assert mu1.reduced_norm() + mu2.reduced_norm() == 2**e*N
    assert θ == mu2 * mu1.conjugate() / N

    #J1 = I*(mu1.conjugate()/N)
    #J2 = I*(mu2.conjugate()/N)
    #assert J1.norm() + J2.norm() == 2**e

    if stats:
        return mu1, mu2, d1, d2, θ, N, I, betaij, num_tries, num_cornacchia
    return mu1, mu2, d1, d2, θ, N, I, betaij, e 

if __name__ == "__main__":
    proof.all(False)
    param = "NIST-I"

    if param == "NIST-I":
        p = 2**248*5 - 1
        e = 246
    elif param == "NIST-III":
        p = 2**376*65 - 1
        e = 374
    elif param == "NIST-V":
        p = 2**500*27 - 1
        e = 498

    assert is_prime(p)
    assert p%4 == 3
    
    B = QuaternionAlgebra(-1, -p)
    #O0 = B.maximal_order()
    # Sage's default is not the isogeny gringos default O0
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i + j) / 2, (1 + k) / 2])
    N = next_prime(randint(p, p**2))

    O_N = O_mod_N(O0, N)

    timings = []
    n_runs = 100

    from tqdm import tqdm
    for _ in tqdm(range(n_runs)):
    #for _ in range(n_runs):
        # print("="*50)
        I = O_N.random_ideal()
        t_start = time()
        qlapoti(I, e, odd_norm_output=True)
        timings.append(time()-t_start)
    print(f"All runs done! Took on average {round(sum(timings)/n_runs, 5)}")
