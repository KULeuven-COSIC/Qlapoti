import sys
sys.path.insert(0, '..')

import random

from sage.all import (
    ZZ,
    ceil, 
    inverse_mod,
    gcd, 
    randint, 
    cached_method,
    Integers,
    Matrix,
    vector,
    factor
)

import time

from norm_eq import qlapoti
from quaternion_helpers.helpers import *
from quaternion_helpers.short_vec import *
from ideal_to_isogeny.dim2_wrapper import Dim2Isogeny
from endomorphisms import *
from theta_structures.couple_point import CouplePoint
from utilities.pairing import weil_pairing_power_two

class IdealToIsogenyQlapotis:
    def __init__(self, J, e_target, precomps, allow_diags = True):
        """
        For a given ideal J, compute the isogeny phi_J : E0 -> EJ
        """

        self.precomps = precomps
        self.representation = None
        self.O0 = precomps.O0
        self.e_target = e_target
        self.J = J
        self.nJ = ZZ(J.norm())
        self._count = (None, None)
        self.allow_diags = allow_diags


    def short_equivalent_ideals(self):
        """
        Given an ideal I of norm N, we find two equivalent ideals I1, I2 of
        coprime odd norms d1, d2 such that such that d1 + d2 = 2 ^ f for some
        f>0
        """
        J = self.J

        e = self.e_target - 2

        mu1, mu2, d1, d2, θ, N, I, betaij, e = qlapoti(J, e, allow_diags = self.allow_diags, odd_norm_output = True)

        assert mu1 in self.O0
        assert mu2 in self.O0
        assert θ in self.O0
        # print(mu1 == make_primitive(mu1, self.O0))
        # print(θ == make_primitive(θ, self.O0))
        self.nI = N
        self.I = I

        return mu1, mu2, d1, d2, e, θ, betaij, e

    def get_representation(self):
        """
        Given an ideal I of norm N, compute the data needed for the
        isogeny. This function only works on the quaternion side.
        """
        if self.representation is None:
            self.representation = self.short_equivalent_ideals()
        return self.representation

    def _precomp_pairings(self, theta, d1, type_end):
        ew0 = self.precomps.eWP0Q0
        P0, Q0 = self.precomps.basis
        
        eP2Q2 = ew0**(theta.reduced_norm())
        eP1Q1 = ew0**(d1**2)

        pair_map = self.precomps.pair_map
        if type_end == 1:
            # Need e(P1, Q1)e(P1, Q2)e(P2, Q1)e(P2, Q2)
            x1, x2, x3, x4 = list(theta)
            eP0Q2 = (
                pair_map['11']**x1 * pair_map['1i']**x2 *
                pair_map['1j']**x3 * pair_map['1k']**x4
            )
            eP2Q0 = (
                pair_map['11']**x1 * pair_map['i1']**x2 *
                pair_map['j1']**x3 * pair_map['k1']**x4
            )

            eP2Q1 = (eP2Q0 ** d1)
            eP1Q2 = (eP0Q2 ** d1)
            w1 = eP1Q1 * eP2Q1
            w2 = eP2Q2 * eP1Q2
            return w1 * w2, [w1, w2]
        
        elif type_end == 2:
            # Need e(iP1, iQ1)*e(iP1, Q2) * e(P2, iQ1) * e(P2, Q2)
            x1, x2, x3, x4 = list(theta)
            eiP0Q2 = (
                pair_map['i1']**x1 * pair_map['ii']**x2 *
                pair_map['ij']**x3 * pair_map['ik']**x4
            )
            eP2iQ0 = (
                pair_map['1i']**x1 * pair_map['ii']**x2 *
                pair_map['ji']**x3 * pair_map['ki']**x4
            )

            eP2iQ1 = (eP2iQ0 ** d1)
            eiP1Q2 = (eiP0Q2 ** d1)
            w1 = eP1Q1 * eP2iQ1
            w2 = eP2Q2 * eiP1Q2
            return w1 * w2, [w1, w2]
        else:
            raise NotImplementedError('aaa')

    def fix_end(self, ker_Phi, f, d1, theta):
        precomps = self.precomps
        P0, Q0 = precomps.basis
        e = precomps.e

        (P1, P2), (Q1, Q2) = ker_Phi

        
        # t0 = time.time()

        prods = []
        f1 = f

        theta2 = tuple(i % 2 for i in list(theta))
        kk = {
            (0, 0, 0, 1):1,
            (0, 1, 1, 1):1,
            (1, 0, 1, 1):2,
            (0, 0, 1, 0):2
        }
        type_end = kk[theta2]
        ends = [type_end]

        if type_end == 1:
            thetap = d1 + theta
            
        else: 
            i,_,_ = self.precomps.B.gens()
            thetap = d1*i + theta
            
        num_prods = thetap.reduced_norm().valuation(2) - 2
        f1 -= 1
        

        
        P1, P2 = evaluate_endomorphism_type(type_end, P1, P2)
        Q1, Q2 = evaluate_endomorphism_type(type_end, Q1, Q2)

        # t1 = time.time()
        # print(f'End time: {t1-t0:.3f}s')


        # Compute pairings (forgetting the multiplication by 2^c)
        _, e_prodtype = self._precomp_pairings(theta, d1, type_end)
        # num_prods = f1 - order_pairing(e_numprods) + 1
        
        # Product steps
        prod_type = None
        if num_prods:
            if e_prodtype[0] ** (2**f1) == e_prodtype[1] ** (2**f1):
                prod_type = 'Q'
            else:
                prod_type = 'P'

        # t2 = time.time()
        # print(f'Prod detect: {t2-t1:.3f}s')
        # print(f"We have to do {num_prods} product(s)...")

        if prod_type == 'P':
            cof2 = 2**(f1 - num_prods)
            K1 = cof2 * P1
            K2 = cof2 * Q2

        elif prod_type == 'Q':
            cof2 = 2**(f1 - num_prods)
            K1 = cof2 * Q1
            K2 = cof2 * P2

        if num_prods > 0:
            K1._order = ZZ(2**num_prods)
            K1._order = ZZ(2**num_prods)
            _phi1 = K1.curve().isogeny(K1)
            _phi2 = K2.curve().isogeny(K2)
            prods.append([_phi1, _phi2])
            P1, Q1 = _phi1(P1), _phi1(Q1)
            P2, Q2 = _phi2(P2), _phi2(Q2)
            f1 -= num_prods

        # t3 = time.time()
        # print(f'Do prod time: {t3 - t2:.3f}s')
        ker_Phi = ((P1, P2), (Q1, Q2))

        # Non-product steps
        Phi = Dim2Isogeny(
                ker_Phi,
                f1,
                extra = self.precomps.extra_torsion_clapotis,
                endomorphisms = ends,
                diagonal_isogenies = prods,
                curves = [P0.curve() for _ in range(2)],
                precomps=self.precomps)

        EI, PI, QI = Phi.embedded_isogeny(
                basis=(P0, Q0),
                degree=d1, e=e, ePQ=precomps.eP0Q0(e))
        # t4 = time.time()
        # print(f'HD time: {t4-t3:.3f}s')
        return EI, PI, QI


    @cached_method
    def images(self):
        """
        Given an ideal I of norm N, this function computes the isogeny
        φI : E0 -> EI.
        Given the 2^e-torsion canonical basis <P0, Q0> we compute φI(P0) and
        φI(Q0)
        """

        precomps = self.precomps
        P0, Q0 = precomps.basis
        e = precomps.e

        # let's find two equivalent ideals of short norm d1 and d2
        # such that d1 + d2 = 2 ^ f
        mu1, mu2, d1, d2, f, theta, betaij, f = self.get_representation()
        
        # let's evaluate θ
        theta_m = precomps.matrix_theta(theta)

        P_, Q_ = self.precomps.evaluate_matrix(theta_m, (P0, Q0))

        # Compute the isogeny corresponding to I1 = I*(mu1_bar)/n(I)
        c = e - f
        P1 = 2 ** c * d1 * P0
        P2 = 2 ** c * P_
        Q1 = 2 ** c * d1 * Q0
        Q2 = 2 ** c * Q_

        ker_Phi = ((P1, P2), (Q1, Q2))
        # Cmon sage
        p = Q1.curve().base().characteristic()

        end_to_fix = all(ci in ZZ for ci in list(theta))

        #assert not end_to_fix

        if end_to_fix:
            t0 = time.time()
            EI, PI, QI = self.fix_end(ker_Phi, f, d1, theta)
            t1 = time.time()
            # print(f'Hd with stuff: {t1 - t0:.3f}s')

        else:
            t0 = time.time()
            Phi = Dim2Isogeny(
                    ker_Phi,
                    f,
                    extra = self.precomps.extra_torsion_clapotis,
                    precomps=self.precomps)

            EI, PI, QI = Phi.embedded_isogeny(
                    basis=(P0, Q0),
                    degree=d1, e=e, ePQ=precomps.eP0Q0(e))
            t1 = time.time()
            # print(f'Only hd: {t1 - t0:.3f}s')

        # TODO: remove
        # ePQ0 = precomps.eP0Q0(e)
        # ePQI = precomps.pairing(PI, QI, e)
        # assert ePQ0 ** d1 == ePQI

        # Now compute J
        theta = mu1 * betaij # Generator of J*I1_bar
        theta /= self.nI

        theta_m = precomps.matrix_theta(theta)
        d1_inv = inverse_mod(d1, 2**e)
        theta_m *= d1_inv

        PJ, QJ = self.precomps.evaluate_matrix(theta_m, (PI, QI))
        # TODO: remove
        # O0 = theta.parent().maximal_order()
        # print(f'{theta in O0 = }')
        # ePQ0 = precomps.eP0Q0(e)
        # ePQJ = precomps.pairing(PJ, QJ, e)
        # assert ePQ0 ** self.nJ == ePQJ

        return EI, PJ, QJ

