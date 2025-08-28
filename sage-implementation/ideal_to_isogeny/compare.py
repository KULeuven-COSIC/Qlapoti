from ideal_to_isogeny_qlapoti import IdealToIsogenyQlapotis
from ideal_to_isogeny_clapotis import IdealToIsogenyClapotis

import sys
sys.path.insert(0, '..')

from applications.PRISM.precomputations_PRSIM import Precomputations_PRISM as Precomputations
from quaternion_helpers.helpers import O_mod_N, reduced_ideal
import time

from sage.all import (
    ZZ, proof, next_prime, randint, QuaternionAlgebra,
    set_random_seed, randint
)

def compare(p, e, c, num_tries):
    rr = randint(1, 2**16)
    # rr = 43746 # 7 products
    print('-'*50)
    print(f'{rr = }')
    print('-'*50)
    set_random_seed(rr)


    precomps = Precomputations(p, c)
    times_qlapoti = []
    times_qlapoti_no_diags = []
    times_clapoti = []

    actual_tries = 0

    while actual_tries < num_tries:

        Isk = precomps.random_ideal()
        Isk = reduced_ideal(Isk)
 
        init_time = time.time()
        phi_I = IdealToIsogenyQlapotis(Isk, e, precomps)
        EI, PI, QI = phi_I.images()
        end_time = time.time()
        times_qlapoti += [end_time - init_time]

        init_time = time.time()
        phi_I = IdealToIsogenyQlapotis(Isk, e, precomps, allow_diags = False)
        EI, PI, QI = phi_I.images()
        end_time = time.time()
        times_qlapoti_no_diags += [end_time - init_time]

        init_time = time.time()
        phi_I = IdealToIsogenyClapotis(Isk, precomps = precomps)
        EI1, PI1, QI1 = phi_I.images()
        end_time = time.time()
        assert EI1.is_isomorphic(EI)
        times_clapoti += [end_time - init_time]
        print(f"Run {actual_tries + 1}")
        actual_tries += 1
        print('=' * 20)

    print("-"*50)
    print(f"Averaged over {num_tries} runs")
    n1 = sum(times_qlapoti)/num_tries
    n2 = sum(times_qlapoti_no_diags)/num_tries
    n3 = sum(times_clapoti)/num_tries
    print(f"IdealToIsogenyQlapotis took {n1:.3f}s")
    print(f"IdealToIsogenyQlapotis (no diags allowed) took {n2:.3f}s")
    print(f"IdealToIsogenyClapotis took {n3:.3f}s")
    print(f"{(n3/n1):.3f}x improvement")
    print(f"{(n3/n2):.3f}x improvement (no diags)")
    print("-" * 50)

if __name__ == "__main__":
    proof.all(False)

    ################################################
    #NIST I
    a = 248
    p = ZZ(5 * 2 ** a - 1)
    c = 122

    # NIST III
    # a = 376
    # p = ZZ(65 * 2**a - 1)
    # c = 122

    # NIST V
    # a = 500
    # p = ZZ(27 * 2**a - 1)
    # c = 122
    num_tries = 100

    compare(p, a, c, num_tries)

