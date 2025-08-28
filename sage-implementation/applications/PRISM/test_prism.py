import time

import params
from sage.all import proof, randint, set_random_seed


from PRISM_sign import PRISM_sign
from precomputations_PRSIM import Precomputations_PRISM

def run_protocol(qlapoti, level=1, precomps_sig=None, precomps_ver=None):
    # r = randint(1, 2**16)
    # set_random_seed(r)
    signer = PRISM_sign('signer', qlapoti=qlapoti, level=level,
                        precomps=precomps_sig)
    verifier = PRISM_sign('verifier', level=level, precomps=precomps_ver)

    init_time = time.time()
    signer.keygen()
    keygen_time = time.time()
    # print(f'Key generation took {(keygen_time - init_time):.3f}s')

    msg = b'Important message ' + str(randint(1, 100)).encode()
    sig = signer.sign(msg)
    sign_time = time.time()
    # print(f'Signing took {(sign_time - keygen_time):.3f}s, counter = {sig[-1]}')

    pk = (signer.pk_curve, signer.pk_sign_basis, signer.eWPpkQpk)
    assert verifier.verify(msg, sig, pk) # Skip for timings

    verif_time = time.time()
    return keygen_time - init_time, sign_time - keygen_time, verif_time - sign_time, verif_time - init_time

if __name__ == '__main__':
    proof.all(False)

    names = ['Keygen', 'Signing', 'Verify', 'Total']
    nruns = 100
    level = 5

    p = params.levels[level]['p']
    precomps_sig = Precomputations_PRISM(p, 122, role = 'signer', a = 'max')
    precomps_ver = Precomputations_PRISM(p, 122, role = 'verifier', a = 'max')


    times_clap = [0 for _ in range(4)]
    times_qlap = [0 for _ in range(4)]
    print(f"Doing {nruns} runs...")
    for k in range(nruns):
        print(f"Run {k+1}")
        r = run_protocol(qlapoti=False, level=level,
                         precomps_sig=precomps_sig, precomps_ver=precomps_ver)
        for i in range(4):
            times_clap[i] += r[i]

    print('Now qlapoti')
    for k in range(nruns):
        print(f"Run {k+1}")
        r = run_protocol(qlapoti=True, level=level)
        for i in range(4):
            times_qlap[i] += r[i]

    print('='*50 + '\nClapoti\n' + '='*50)
    for i in range(4):
        tt = times_clap[i] / nruns
        print(f'{names[i]} time: {tt:.3f}s')
    print('='*50 + '\nQlapoti\n' + '='*50)
    for i in range(4):
        tt = times_qlap[i] / nruns
        print(f'{names[i]} time: {tt:.3f}s')
    print('='*50)
    n2 = times_clap[0] / nruns
    n1 = times_qlap[0] / nruns
    print(f"Improvement to KeyGen {(n2/n1):.3f}x !!")
    n2 = times_clap[1] / nruns
    n1 = times_qlap[1] / nruns
    print(f"Improvement to Signing {(n2/n1):.3f}x !!")
    print('='*50)
    



