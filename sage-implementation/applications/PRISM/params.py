################################################
# parameters for the PRISM-{id,sig} scheme
# Semi-Mersenne prime for NIST I
a = 248 # Available 2-torsion
p = 5* 2**a - 1

# Respective size of the responses
id_torsion = 136 # lambda + log(lambda)
sign_torsion = 219
sign_torsion_sidh = 247

# Compress signature (always True for the id scheme)
point_compression_sign = False

levels = {
    1 : {
        'p':5*2**248 - 1,
        'a':248,
        'sign_torsion':219,
        'id_torsion':136
    },
    3: {
        'p':65*2**376 - 1,
        'a':376,
        'sign_torsion':376,
        'id_torsion':200
    },
    5: {
        'p':27*2**500 - 1,
        'a':500,
        'sign_torsion':500,
        'id_torsion':265
    }
}


################################################
# parameters for the PRISM-{id,sig} scheme using SIDH

# primes of SIDH shape
p1 = 13 * 2**126 * 3**78  - 1 
p2 = 5  * 2**193 * 3**122 - 1 
p3 = 11 * 2**257 * 3**163 - 1



################################################

# Parameters for the push even scheme
param_push_even_lvl1 = {
    'p': 65 * 2**376- 1,
    'e': 376, # full torsion
    'a': 248, # sign torsion
    'u': 128, # push torsion
    'hash_len': 248,
    'num_steps': 2,
}

param_push_even = param_push_even_lvl1

################################################

