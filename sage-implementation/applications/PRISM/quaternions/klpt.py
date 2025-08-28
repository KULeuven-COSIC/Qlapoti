"""

WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

THIS FILE WILL NOT WORK PROPERLY AS THE OLD KLPT LIKE CODE HAD A BUNCH OF
GLOBALS BECAUSE PAST JACK MADE VERY BAD CODING DECISIONS BUT IT SHOULD BE ABLE
TO B FIXED WITH SMALL CHANGES

WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

"""


# Python imports
from random import choice

# Sage imports
from sage.all import gcd, ZZ, floor, sqrt, randint, kronecker, random_prime

# Local imports
from ideals import (
    chi,
    reduced_basis,
    quadratic_norm,
)
from utilities.sum_of_squares import sum_of_squares
from lattices import generate_small_norm_quat


def is_inert(ω, p):
    """
    Given an element ω ∈ B, and a prime p, returns if
    p is inert in ℤ[ω]

    For ℤ[ω], a prime is inert when the Kronecker symbol:
        (-d / p) = -1
    """

    def discriminant(ω):
        """
        Computes the discriminant of an element of a
        quaternion algebra, denoted Δ(ω)
        """
        Trd = ω.reduced_trace()
        Nrd = ω.reduced_norm()
        if Trd not in ZZ or Nrd not in ZZ:
            raise ValueError(f"The element ω must be integral")

        return ZZ(Trd**2 - 4 * Nrd)

    def recover_d(ω):
        """
        Given an integral element of a Quaternion Algebra,
        and its discriminant Δ(ω), computes the integer d
        such that:

        Δ(ω) = -d    if d ≡ 3 mod 4
             = -4d   otherwise
        """
        Δ = discriminant(ω)
        if Δ % 4 == 0:
            return -(Δ // 4)
        return -Δ

    d = recover_d(ω)
    return kronecker(-d, p) == -1


def inert_prime(bound, d):
    """
    Input: An upper bound for the output prime
           d, the integer such that ω^2 = -d for
           ℤ[ω]

    Output: A prime < bound which is inert in
            the ring ℤ[ω]
    """
    while True:
        p = random_prime(bound)
        if kronecker(-d, p) == -1:
            return p


def generate_small_norm_quat_random(Ibasis, coeff_bound, search_bound):
    """
    Pick a random linear combination from Ibasis to compute an element
    α ∈ B_{0, ∞}
    """
    for _ in range(search_bound):
        xs = [randint(-coeff_bound, coeff_bound) for _ in range(len(Ibasis))]
        if gcd(xs) != 1:
            continue
        α = sum(αi * x for αi, x in zip(Ibasis, xs))

        yield α


# TODO: prime_norm_heuristic should be derived, before it was a global
def prime_norm_algebra_element(
    nI,
    Ibasis,
    coeff_bound,
    search_bound,
    previous=set(),
    allowed_factors=None,
    random_elements=False,
    prime_norm_heuristic=200,  # TODO this should be set as we do in SQISign-SageMath
):
    """
    Find an element α ∈ B_{0, ∞} with small,
    prime scaled norm.

    Optional: `allowed_factors` allows the norm to
               be composite, where it is expected that
               the result is a large prime multiplied by
               small factors dividing allowed_factors.
    """

    if random_elements:
        small_elements = generate_small_norm_quat_random(
            Ibasis, coeff_bound, search_bound
        )
    else:
        max_norm_bound = prime_norm_heuristic
        small_elements = generate_small_norm_quat(
            Ibasis, max_norm_bound, count=search_bound
        )

    for α in small_elements:
        α_norm = ZZ(α.reduced_norm()) // nI

        # Even norms can be rejected early
        # as allowed_factors is either odd
        # or None.
        if α_norm % 2 == 0:
            continue

        # We can allow α to have composite norm
        # if small factors are within the KLPT
        # target norm T.
        α_norm_reduced = α_norm
        if allowed_factors:
            g = gcd(α_norm, allowed_factors)
            if g != 1:
                α_norm_reduced //= g

        # If we've failed with this norm before
        # continue
        if α_norm_reduced in previous:
            continue

        # Check if the element has prime norm
        if α_norm_reduced.is_pseudoprime():
            # Check if the prime small enough
            if α_norm_reduced < prime_norm_heuristic:
                return α, α_norm
    return None, None


# TODO: logp should be derived, before it was a global
def EquivalentPrimeIdealHeuristic(
    I, previous=set(), allowed_factors=None, random_elements=False, logp=256
):
    """
    Given an ideal I with norm nI, attempts
    to find an equivalent ideal J with prime norm.

    If unsuccessful, returns None
    """
    # TODO: what's a good initial small search?
    coeff_bound = max((floor(logp / 10)), 7)

    # TODO: what's a good search bound?
    search_bound = max(coeff_bound**4, 4096)

    # Norm of Ideal
    nI = ZZ(I.norm())

    # Compute the Minkowski reduced basis
    Ibasis = reduced_basis(I)

    # Find an element with small prime norm
    α, N = prime_norm_algebra_element(
        nI,
        Ibasis,
        coeff_bound,
        search_bound,
        previous=previous,
        allowed_factors=allowed_factors,
        random_elements=random_elements,
    )

    if α is None:
        print(f"DEBUG [EquivalentPrimeIdealHeuristic] No equivalent prime found")
        return None, None, None

    assert ZZ(α.reduced_norm()) // nI == N
    assert α in I

    # Compute the ideal given α
    J = chi(α, I)

    return J, N, α


def RepresentIntegerHeuristic(M, parity=False):
    """
    Algorithm 1 (Page 8)

    Given an integer M, with M > p, attempts to
    find a random element γ with norm M.

    If no element is found after `bound` tries,
    returns none
    """

    def RepresentInteger(M, z, t, parity=False):
        M_prime = M - p * quadratic_norm(z, t)
        two_squares = sum_of_squares(M_prime)
        if two_squares:
            x, y = two_squares
            if parity and (x + t) % 2 == 0 and (y + z) % 2 == 0:
                return None
            return x + ω * y + j * (z + ω * t)
        # No solution for the given M
        return None

    if M <= p:
        raise ValueError(f"Can only represent integers M > p.")
    m = max(floor(sqrt(M / (p * (1 + q)))), 5)

    # TODO: how many times should we try?
    for _ in range(m**2):
        z = randint(-m, m)
        t = randint(-m, m)
        γ = RepresentInteger(M, z, t, parity=parity)

        if γ is not None:
            # Found a valid solution, return
            assert γ.reduced_norm() == M, "The norm is incorrect"
            assert γ in O0, "The element is not contained in O0"
            return γ

    # No solution found, return None
    print(f"DEBUG [RepresentIntegerHeuristic]: No solution found")
    return None


def EquivalentRandomEichlerIdeal(I, Nτ):
    """
    Algorithm 6 (SQISign paper)

    Input:  I a left O-ideal
    Output: K ∼ I of norm coprime with Nτ
    """
    nI = I.norm()
    # TODO: what should the size of `bound` be
    bound = 10

    # Step 1: find an element ωS such that Nτ is inert in ℤ[ωS]
    O = I.left_order()
    while True:
        ωS = sum([randint(-bound, bound) * b for b in O.basis()])
        if is_inert(ωS, Nτ):
            break

    # Step 2: find a random element γ in I such that n(γ)/n(I) is coprime with Nτ
    while True:
        γ = sum([randint(-bound, bound) * b for b in I.basis()])
        if gcd(γ.reduced_norm() // nI, Nτ) == 1:
            break

    # Step 3: select a random class (C : D) ∈ P1(Z/Nτ Z).
    x = randint(0, Nτ)
    if x == p:
        C, D = 1, 0
    else:
        C, D = x, 1

    # Step 4: set β = (C + ωSD)γ.
    β = (C + ωS * D) * γ

    # Step 5: return K = χI(β)
    return chi(β, I)
