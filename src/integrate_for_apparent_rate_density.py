'''
We wrote a general equation for the apparent rate density, as a function of the
apparent planet radius r_a and apparent mass M_a.

Here, we evaluate it for particular "true" (volume-limited) occurrence rate
densities.

The models are as follows:
/Model #1: twin binary, same planet

XModel #2: varying q, same planet (WARNING: not implemented)

/Model #3: varying q, powerlaw+const r
/Model #4: varying q, r with a gap
/Model #5: twin binary, powerlaw r
/Model #6: varying q, powerlaw r

Example usage, for each model type:
>>> python integrate_for_apparent_rate_density.py --modelnumber 1 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 1
>>> python integrate_for_apparent_rate_density.py --modelnumber 2 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python integrate_for_apparent_rate_density.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python integrate_for_apparent_rate_density.py --modelnumber 4 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python integrate_for_apparent_rate_density.py --modelnumber 5 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 22.5
>>> python integrate_for_apparent_rate_density.py --modelnumber 6 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5

'''
from __future__ import division, print_function

import numpy as np, pandas as pd
from scipy import signal
from scipy.integrate import trapz
import argparse

global α, β, γ, δ
α = 3.5
β = 0
γ = 0
δ = -2.92
global Z_0, Z_1, Z_2
Z_0, Z_1 = 0.5, 0.5
global BF, model_number, r_pu, r_pl
r_pl = 2 # [Rearth]. Lower bound for truncation.


def Gamma_0(r, M, Z_0, model_number, norm_r=None):
    '''
    occ rate density around singles
    '''

    if model_number in [1,2]:
        if r == 1:
            Γ_0 = Z_0
        else:
            Γ_0 = 0

    elif model_number == 3:
        f_r = min( r**δ, r_pl**δ )
        Γ_0 = Z_0 * f_r / norm_r

    elif model_number in [5,6]:
        Γ_0 = Z_0 * r**δ / norm_r

    elif model_number == 4:
        f_r = 0 if (r > 1.5) & (r < 2) else min( r**δ, r_pl**δ )
        Γ_0 = Z_0 * f_r / norm_r

    return Γ_0


def Gamma_1(r, M, Z_1, model_number, prefactor=None, norm_r=None, q_grid=None):
    '''
    occ rate density around primaries

    prefactor:
        e.g., "1/A = sqrt(2)" in Model #1 (needed for delta function norm)
    '''

    if model_number in [1]:

        r_0 = 1
        # janky delta function
        Γ_1 = Z_1/prefactor * (np.isclose(r,r_0,atol=1e-10)).astype(int)

    elif model_number in [2]:

        raise NotImplementedError

        r_0 = 1
        # δ(g(q)) = sum over roots of δ(q-q_i) / |g'(q_i)|
        g_of_q = r - r_0
        for i, eps in enumerate(np.diff(q_grid)):
            g_prime_of_q = (g_of_q[i+1] - g_of_q[i])/eps
        #TODO: I currently don't know how to take this numerical derivative
        #and/or root-find in a reliable way for this.
        _ = np.argmin( abs(g_of_q) )
        Γ_1 = Z_1/prefactor * (np.isclose(r,r_0,atol=1e-10)).astype(int)

    elif model_number == 3:
        f_r = np.minimum( r**δ, r_pl**δ )
        Γ_1 = Z_1 * f_r / norm_r

    elif model_number in [5,6]:
        Γ_1 = Z_1 * r**δ / norm_r

    elif model_number == 4:
        f_r = np.minimum( r**δ, r_pl**δ )
        f_r[(r > 1.5) & (r < 2)] = 0
        Γ_1 = Z_1 * f_r / norm_r

    return Γ_1


def Gamma_2(r, M, Z_2, model_number, prefactor=None, norm_r=None):
    '''
    occ rate density around secondaries
    '''

    if model_number in [1,2]:
        r_0 = 1
        Γ_2 = Z_2/prefactor * (np.isclose(r,r_0,atol=1e-10)).astype(int)

    elif model_number == 2:
        raise NotImplementedError
        #same story as for Gamma_1, except now the multiple roots matters.

    elif model_number == 3:
        f_r = np.minimum( r**δ, r_pl**δ )
        Γ_2 = Z_2 * f_r / norm_r

    elif model_number in [5,6]:
        Γ_2 = Z_2 * r**δ / norm_r

    elif model_number == 4:
        f_r = np.minimum( r**δ, r_pl**δ )
        f_r[(r > 1.5) & (r < 2)] = 0
        Γ_2 = Z_2 * f_r / norm_r

    return Γ_2


def A(q):
    return 1/np.sqrt(1+q**α)


def B(q):
    return (1/q) * 1/np.sqrt(1 + q**(-α))


def f_of_q(q, model_number):
    '''
    volume-limited binary mass ratio distribution.
    '''

    if model_number in [1,5]:
        f_q = signal.unit_impulse(len(q), -1)

    elif model_number in [2,3,4,6]:
        f_q = q**β

    norm_q = trapz(f_q, q)
    f_q /= norm_q
    return f_q, norm_q


def run_unit_tests(
        model_number, Γ_a, r_a, Γ_0,
        Z_0, Z_1, Z_2, BF, r_pu
        ):

    if model_number == 1:

        μ = BF/(1-BF) * 2**(3/2)

        r_p = 1

        assert Γ_a[r_a==r_p] == Z_0 / (1+μ)

        assert np.isclose( float(Γ_a[np.isclose(r_a, r_p*2**(-1/2))]),
                           (Z_1+Z_2) * μ / (1+μ),
                           rtol=1e-12)

    elif model_number == 2:

        print('analytic comparison not yet working')
        import IPython; IPython.embed()
        raise NotImplementedError

        q_grid = np.arange(1e-6, 1+1e-6, 1e-6)
        prob_vl_q = q_grid**β
        norm_q = trapz(prob_vl_q, q_grid)

        # Derived 2017/11/30.1. Has the same functional form as the one derived
        # on 2017/10/10.1 (but has explicit normalization).
        r_p = 1
        I_1 = Z_1/norm_q * 2/α * r_p/(r_a**2) * ( (r_p/r_a)**2 -1 )**((β-α+1)/α) * (r_p/r_a)**4
        I_1[r_a == 0] = 0
        I_1[r_a == 1] = 0

        μ = get_mu(BF, model_number)
        Γ_primaries_analytic = 1/(1+μ) * BF/(1-BF) * I_1

        # Expect everything from r_p/sqrt(2) to just below r_p to be from the
        # primaries. Going down to 0.73, not 1/sqrt(2) because the funky r_a(q)
        # relation for the i=2 case goes ABOVE sqrt(2) (!).
        vals_num = Γ_a[(r_a>0.73) & (r_a<0.99)]
        vals_anal = Γ_primaries_analytic[(r_a>0.73) & (r_a<0.99)]

        # Assert the numeric and analytic distributions agree to 0.5%.
        assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=5e-3)


    elif model_number == 3:

        ##################################################################
        # check ratio of inferred and single rate densities at r_a>r_pl. #
        ##################################################################
        q = get_q_grid(model_number)
        f_q = q**β
        norm_q = trapz(f_q, q)
        X_num = Γ_a/Γ_0
        μ = get_mu(BF, model_number)
        if Z_0 == Z_1 == Z_2:
            assert np.isclose(
                    np.mean(X_num[(r_a>r_pl) & (r_a<=r_pu)]),
                    1/(1+μ) * (1 + 1/norm_q * BF/(1-BF) * (1.41425)),
                    rtol=1e-6
                    )

        ##################################################################
        # check ratio of inferred and single rate densities at r_a<r_pl. #
        ##################################################################
        _grid = np.arange(1e-3, r_pu+1e-3, 1e-3)
        prob_r = np.minimum( _grid**δ, r_pl**δ )
        norm_r = trapz(prob_r, _grid)

        sel_inds = r_a<r_pl
        r_a = r_a[sel_inds]

        # First, get I_1.

        # minimization for dummies:
        # get q_l_of_ra precise to 2e-6 (or whatever the q_grid spacing is) by
        # just taking the minimum of the absolute value.
        q_l_of_ra = q[np.argmin(np.abs(
                      ((1 + q**α))[:,None] - ((r_pl/r_a)**2)[None,:]
                      ),
                      axis=0)]

        # there's definitely a way to vectorize this
        int1_1, int1_2 = [], []
        for this_q_l_of_ra, this_ra in zip(q_l_of_ra, r_a):
            int1_1.append(trapz(
                    q[q>this_q_l_of_ra]**β * (1 + q[q>this_q_l_of_ra]**α)**((δ+4)/2),
                    q[q>this_q_l_of_ra]
                    ))
            int1_2.append(trapz(
                    q[q<this_q_l_of_ra]**β * (1 + q[q<this_q_l_of_ra]**α)**(2),
                    q[q<this_q_l_of_ra]
                    ))
        int1_1, int1_2 = np.array(int1_1), np.array(int1_2)

        I_1 = Z_1 / (norm_q * norm_r) * (
              r_a**δ * (int1_1) +
              r_pl**δ * (int1_2)
              )

        # Second, get I_2.
        q_t = 0.9210924370 # turning point
        minval = 1.9796263301 # q^2 * (1 + q^{-\alpha}) at q_t
        d_inds = (minval < (r_pl/r_a)**2) & ((r_pl/r_a)**2 < 2) # double-valued

        I_2 = np.zeros_like(I_1)

        # subcase #1, when NOT double-valued
        q_l1_of_ra = q[np.argmin(np.abs(
                        (q**2 * (1 + q**(-α)))[:,None] - (((r_pl/r_a)**2)[~d_inds])[None,:] ),
                       axis=0)]

        # there's definitely a way to vectorize this
        int2_1, int2_2 = [], []
        for this_q_l1_of_ra, this_ra in zip(q_l1_of_ra, r_a[~d_inds]):

            _q = q[q<this_q_l1_of_ra]

            int2_1.append(
                    trapz(
                    _q**(β+γ+δ+8/3) * (1 + _q**α)**(3/2) * (1 + _q**(-α))**((δ+1)/2),
                    _q
                    ))

            _q = q[q>this_q_l1_of_ra]

            int2_2.append(
                    trapz(
                    _q**(β+γ+8/3) * (1 + _q**α)**(3/2) * (1 + _q**(-α))**(1/2),
                    _q
                    ))

        int2_1, int2_2 = np.array(int2_1), np.array(int2_2)

        I_2[~d_inds] += Z_1 / (norm_q * norm_r) * (
                        r_a[~d_inds]**δ * (int2_1[~d_inds]) +
                        r_pl**δ * (int2_2[~d_inds])
                        )

        # subcase #2, when double-valued
        if np.any(d_inds):
            raise NotImplementedError
            I_2[d_inds] += 0 # todo

        # Finally, use them compute Γ_a, apparent rate density
        Γ_a_from_math = (1/(1+μ))*(
                Γ_0[sel_inds] + BF/(1-BF)*(I_1 + I_2)
              )

        if Z_0 == Z_1 == Z_2:
            assert np.all(np.isclose(
                Γ_a[sel_inds][r_a < r_pl*2**(-1/2)],
                Γ_a_from_math[r_a < r_pl*2**(-1/2)],
                rtol=5e-6))


    elif model_number == 5:

        X_num = Γ_a/Γ_0
        μ = get_mu(BF, model_number)
        assert np.isclose(
                np.mean(X_num[(r_a>2) & (r_a<=r_pu/2**(1/2))]),
                (1+2**((δ+3)/2)*μ)/(1+μ),
                rtol=1e-12
                )


    elif model_number == 6:

        # check single rate density agrees with analytic value
        r = r_a
        norm_r = trapz(r**δ, r)
        Γ_vl = r**δ * Z_0 / norm_r
        np.testing.assert_allclose(
                Γ_vl/Γ_0, np.ones_like(Γ_vl/Γ_0),
                rtol=1e-12)

        # check ratio of inferred and single rate densities
        q = get_q_grid(model_number)
        f_q = q**β
        norm_q = trapz(f_q, q)
        X_num = Γ_a/Γ_0
        μ = get_mu(BF, model_number)
        assert np.isclose(
                np.mean(X_num[(r_a>1) & (r_a<=r_pu/2**(1/2))]),
                1/(1+μ) * (1 + 1/norm_q * BF/(1-BF) * (1.41425)),
                rtol=1e-6
                )


def get_mu(BF, model_number):

    if model_number in [1,5]:
        μ = BF/(1-BF) * 2**(3/2)

    elif model_number in [2,3,4,6]:
        I = 0.4841741 # see derivation pdf
        μ = BF/(1-BF) * ( 2**(3/2) - 3*I )

    return μ


def Gamma_a(r_a, M_a, f_q, q, A_q, B_q, μ, norm_r=None, norm_q=None,
        debugging=False):

    Γ_0 = Gamma_0(r_a, M_a, Z_0, model_number, norm_r=norm_r)

    #NOTE: for model #2, this prefactor is wrong. need to implement
    #derivative? or something else to capture the noramlization funkiness.
    Γ_1 = Gamma_1(
            r_a/A_q, M_a,
            Z_1,
            model_number,
            prefactor=1/A_q,
            norm_r=norm_r,
            q_grid=q)

    I_1 = trapz(f_q * (A_q**(-4)) * Γ_1, q, axis=0)

    Γ_2 = Gamma_2(
            r_a/B_q, q*M_a,
            Z_2,
            model_number,
            prefactor=1/B_q,
            norm_r=norm_r)

    I_2 = trapz((q**(5/3)) * \
                f_q * \
                (A_q**(-3)) * \
                (B_q**(-1)) * \
                Γ_2, \
                q, axis=0)

    Γ_a = (1/(1+μ))*(
            Γ_0 + BF/(1-BF)*(I_1 + I_2)
          )

    if debugging:

        if model_number == 6:
            # int_0^1 (dq q**β * (1+q**α)**((δ+4)/2)) = 1.10812
            # (171211_model3_integrals.nb)
            I_1_from_math = Z_1/(norm_r*norm_q) * r_a**δ * 1.10812
            I_2_from_math = Z_2/(norm_r*norm_q) * r_a**δ * 0.306134

            print('{:.4g}\n{:.4g}'.format(
                (I_1-I_1_from_math)/I_1, (I_2-I_2_from_math)/I_2)
                )

            Γ_a_from_math = (1/(1+μ))*(
                    Γ_0 + BF/(1-BF)*(I_1_from_math + I_2_from_math)
                  )
            print('{:.4g}'.format(
                (Γ_a-Γ_a_from_math)/Γ_a)
                )

            #NOTE: every point should match. Your numerical integration agrees with
            #your analytic integration for Γ_a. Ditto for Γ_0.
            #import IPython; IPython.embed()
            Γ_0_from_math = r_a**δ * Z_0 / (norm_r)
            X_from_math = Γ_a_from_math/Γ_0_from_math
            X_num = Γ_a/Γ_0

        elif model_number == 2 and r_a>0.73 and r_a<0.99:
            # int_0^1 (dq q**β * (1+q**α)**((δ+4)/2)) = 1.10812
            # (171211_model3_integrals.nb)
            r_p = 1
            I_1_from_math = Z_1/norm_q * 2/α * r_p/(r_a**2) * ( (r_p/r_a)**2 -1 )**((β-α+1)/α) * (r_p/r_a)**4

            print('{:.4g}'.format(
                (I_1-I_1_from_math)/I_1)
                )

            #NOTE: every point should match. Your numerical integration agrees with
            #your analytic integration for Γ_a. Ditto for Γ_0.
            import IPython; IPython.embed()

            Γ_a_from_math_ish = (1/(1+μ))*(
                    Γ_0 + BF/(1-BF)*(I_1_from_math + I_2)
                  )
            print('{:.4g}'.format(
                (Γ_a-Γ_a_from_math_ish)/Γ_a)
                )

            Γ_0_from_math = r_a**δ * Z_0 / (norm_r)
            X_from_math = Γ_a_from_math/Γ_0_from_math
            X_num = Γ_a/Γ_0

    return Γ_a, Γ_0


def get_apparent_radius_grid(model_number):

    if model_number == 1:
        r_a_grid = np.array([2**(-1/2),1])

    elif model_number == 2:
        r_a_grid = np.arange(5e-3,1+5e-3,5e-3)

    elif model_number in [5]:
        r_a_grid = np.arange(1e-2,r_pu+1e-2,1e-2)

    elif model_number in [3, 4, 6]:
        r_a_grid = np.arange(1e-1,r_pu+1e-1,1e-1)

    return r_a_grid


def get_q_grid(model_number):

    if model_number in [1,5]:
        q = np.arange(0.1,1+0.1,0.1)

    elif model_number in [2,3,4,6]:
        # 1e-5 needed for any non-resolution limited result
        q = np.arange(2e-6,1+2e-6,2e-6)

    return q


def get_apparent_rate_density(r_a_grid, M_a_grid, model_number):
    # gives the rate density around singles too

    μ = get_mu(BF, model_number)

    q = get_q_grid(model_number)

    f_q, norm_q = f_of_q(q, model_number)

    if model_number in [5,6]:
        norm_r = trapz(r_a_grid**δ, r_a_grid)
    elif model_number == 3:
        r_grid = np.arange(1e-3, r_pu+1e-3, 1e-3)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        norm_r = trapz(prob_r, r_grid)
    elif model_number == 4:
        r_grid = np.arange(1e-3, r_pu+1e-3, 1e-3)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        prob_r[(r_grid > 1.5) & (r_grid < 2)] = 0
        norm_r = trapz(prob_r, r_grid)
    elif model_number == 1:
        norm_r = 1

    A_q = A(q)
    B_q = B(q)

    Γ_a, Γ_0 = [], []
    print('ind/max_ind, rvalue')
    for r_a_ind, r_a in enumerate(r_a_grid):
        for M_a in M_a_grid:
            if r_a_ind % max(1,len(r_a_grid)//10) == 0:
                print('{:d}/{:d}, {:.2f}'.format(r_a_ind,len(r_a_grid),r_a))
            _0, _1 = Gamma_a(r_a, M_a, f_q, q, A_q, B_q, μ, norm_r=norm_r,
                    norm_q=norm_q)
            Γ_a.append(_0)
            Γ_0.append(_1)

    Γ_a = np.array(Γ_a)
    Γ_a = Γ_a.reshape(len(r_a_grid),len(M_a_grid)).flatten()
    Γ_0 = np.array(Γ_0)
    Γ_0 = Γ_0.reshape(len(r_a_grid),len(M_a_grid)).flatten()

    return Γ_a, Γ_0


def write_output(df):

    savdir = '../data/numerics/'
    outname = 'integrate_model_{:d}_Zsub2_{:.2f}_rpu_{:.2f}.out'.format(
            model_number, Z_2, r_pu
            )
    df.to_csv(savdir+outname, index=False)
    print('saved to {:s}'.format(savdir+outname))



if __name__ == '__main__':

    ##########################################
    parser = argparse.ArgumentParser(
        description='Integrate the general equation for apparent rate density.')

    parser.add_argument('-mn', '--modelnumber', type=int, default=None,
        help='1, 2 or 3.')
    parser.add_argument('-ztwo', '--ZsubTwo', type=float, default=None,
        help='integrated occ rate for secondaries')
    parser.add_argument('-bf', '--binaryfrac', type=float, default=None,
        help='BF = n_b/(n_s+n_b), for n number density')
    parser.add_argument('-rpu', '--upperradiuscutoff', type=float, default=None,
        help='the maximum allowed planet radius in units of Rearth')

    args = parser.parse_args()

    model_number = args.modelnumber
    Z_2 = args.ZsubTwo
    BF = args.binaryfrac
    r_pu = args.upperradiuscutoff
    ##########################################

    M_a_grid = np.array([1.0])

    r_a_grid = get_apparent_radius_grid(model_number)

    Γ_a, Γ_0 = get_apparent_rate_density(r_a_grid, M_a_grid, model_number)

    df = pd.DataFrame({'r':r_a_grid, 'Γ_a':Γ_a, 'Γ_0':Γ_0})

    write_output(df)

    run_unit_tests(
            model_number, Γ_a, r_a_grid, Γ_0,
            Z_0, Z_1, Z_2, BF, r_pu
            )
