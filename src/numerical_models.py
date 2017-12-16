# -*- coding: utf-8 -*-
'''
See README.md for full description.

The models are as follows:
Model #1: twin binary, same planet
Model #2: varying q, same planet
Model #3: varying q, varying r
Model #4: varying q, r with a gap
Model #5: twin binary, varying r

Example usage, for each model type:
>>> python numerical_models.py --modelnumber 1 --ZsubTwo 0.5 --binaryfrac 0.1
>>> python numerical_models.py --modelnumber 2 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python numerical_models.py --modelnumber 4 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
>>> python numerical_models.py --modelnumber 5 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 22.5

Alternatively, use a wrapper like `run_numerical_models.py`.
'''

from __future__ import division, print_function

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from math import pi as π
import os, argparse, decimal


def check_inputs(γ, Z_2, model_number, BF, r_pu):
    '''
    Verify that input values are reasonable, otherwise raise errors.
    '''

    if γ != 0:
        raise NotImplementedError

    if abs(decimal.Decimal(str(Z_2)).as_tuple().exponent) > 2:
        # File names assume 2 digits of precision for Z_2.
        raise NotImplementedError

    if model_number > 5 or model_number < 1:
        raise NotImplementedError

    if not (BF > 0):
        raise AssertionError

    if model_number in [2,3,4]:
        if not r_pu > 0:
            raise AssertionError


def make_stellar_population(quickrun, model_number, BF, α, β, γ, δ):
    '''
    STELLAR POPULATION

    Construct a pandas DataFrame in which each row is a star that has been
    "selected". Each selected star has properties:

        * star_type (str): single, primary, or secondary

        * is_searchable (bool): whether the star is searchable (determined by
          completeness for star-specific q).

        * q (float): binary mass ratio, if applicable

    Ensure:

        * primary and secondary in same system get same q
        * 0 < q < 1
        * there are as many selected primaries and secondaries

    RETURNS:

        N_0, N_1, N_2: number of selected singles, primaries, secondaries
        is_searchable: boolean array, indices matched to df, of searchability
                       per star
        df: pandas dataframe with all stellar data.
        q: mass ratios for binaries (length N_1).
    '''

    # arbitrary number of selected single stars. 1e6 is pretty fast.
    N_0 = int(5e3) if quickrun else int(5e6)

    B = BF / (1-BF) # prefactor in definition of μ

    # Get number of selected primaries.
    if model_number == 1 or model_number == 5:
        q = 1
        μ = B * (1 + q**α)**(3/2)
        N_1 = int(N_0 * μ) # number of selected primaries
    elif model_number == 2 or model_number == 3 or model_number == 4:
        integral = 0.484174086070513 # 17/10/14.2 analytic result
        N_d = N_0 * B * (2**(3/2) - 3*integral)
        N_1 = int(N_d) # number of primaries = number of double star systems

    N_2 = N_1 # number of selected secondaries

    # Construct population of stars: singles, primaries, secondaries.
    df = pd.concat((
        pd.DataFrame({'star_type':[type_i for _ in range(N_i)]}) for N_i,type_i in
        [(N_0,'single'),(N_1,'primary'),(N_2,'secondary')]
        ), ignore_index=True )

    # Assign values of `q` to primaries and secondaries
    if model_number == 1 or model_number == 5:
        df['q'] = q
        df.loc[df['star_type'] == 'single', 'q'] = np.nan
        q = np.ones(N_1)

    elif model_number == 2 or model_number == 3 or model_number == 4:
        df['q'] = np.nan

        # Inverse transform sampling to get samples of q
        q_grid = np.arange(1e-6, 1+1e-6, 1e-6)
        prob_ml_q = q_grid**β * (1+q_grid**α)**(3/2)
        prob_ml_q /= trapz(prob_ml_q, q_grid)
        cdf_ml_q = np.append(0, np.cumsum(prob_ml_q)/np.max(np.cumsum(prob_ml_q)))
        func = interp1d(cdf_ml_q, np.append(0, q_grid))
        q_samples = func(np.random.uniform(size=int(N_1)))

        df.loc[df['star_type'] == 'primary', 'q'] = q_samples
        df.loc[df['star_type'] == 'secondary', 'q'] = q_samples

        q = q_samples

    np.testing.assert_equal(
            len(df[df['star_type']=='primary']),
            len(df[df['star_type']=='secondary'])
            )

    np.testing.assert_array_equal(
            df[df['star_type']=='primary']['q'],
            df[df['star_type']=='secondary']['q'])

    np.testing.assert_array_less(
            np.zeros(len(q)),
            q
            )

    # Select stars that are searchable. This means they win the completeness
    # lottery. NB this is independent of geometric transit probability.
    Q_c0 = 1
    single_is_searchable = np.ones(N_0).astype(bool)

    # NOTE: the inclusion of this "searchability" draw is outdated. (In drafts
    # before mid-November 2017, I made an error in counting searchable
    # binaries -- I thought that e.g., 1 in 8 selected twin binaries would be
    # searchable. However, a star is selected based on its flux and the
    # planet's *apparent* radius, rather than the planet's true radius. Thus
    # all selected twin binaries are searchable. (Per KM's 2017/11/21 note).
    Q_c1 = 1 #(1+q**α)**(-3)
    primary_is_searchable = ( np.random.rand((N_1)) < Q_c1 )

    Q_c2 = 1 #(1+q**(-α))**(-3) * q**(-5)
    secondary_is_searchable = ( np.random.rand((N_2)) < Q_c2 )

    is_searchable = np.concatenate(
        (single_is_searchable,  primary_is_searchable, secondary_is_searchable)
        )
    df['is_searchable'] = is_searchable

    _1 = len(df[(df['is_searchable']==True)&(df['star_type']=='primary')])/\
            len(df[df['star_type']=='primary'])
    print("population's completeness fraction for primaries: {:.3g}".
            format(_1))

    _2 = len(df[(df['is_searchable']==True)&(df['star_type']=='secondary')])/\
            len(df[df['star_type']=='secondary'])
    print("population's completeness fraction for secondaries: {:.3g}".
            format(_2))

    return N_0, N_1, N_2, is_searchable, df, q


def make_planet_population(model_number, df, q, is_searchable, Z_0, Z_1, Z_2,
        N_0, N_1, N_2, δ, r_pu, debugging=False):
    '''
    Given the dataframe of stars, assign planets. Further, assign whether they
    transit, and whether they are detected.

    kwargs:
        debugging: True/False if debugging/not debugging.
    '''

    #######################################################
    # Select stars with planets, and assign planet radii. #
    #######################################################
    single_has_planet = ( np.random.rand((N_0)) < Z_0 )
    primary_has_planet = ( np.random.rand((N_1)) < Z_1 )
    secondary_has_planet = ( np.random.rand((N_2)) < Z_2 )

    has_planet = np.concatenate(
        (single_has_planet,  primary_has_planet, secondary_has_planet)
        )
    df['has_planet'] = has_planet

    if model_number == 1 or model_number == 2:
        # Assign arbitrary planet radius.
        r_p = 1
        r_pu = r_p
        df['r'] = r_p
        df.loc[df['has_planet'] == False, 'r'] = np.nan

    elif model_number == 3 or model_number == 5:
        r_pl = 2 # [Rearth]. Lower bound for truncation.

        # Inverse transform sample to get radii. Drawing from powerlaw
        # distribution above r_pl, and constant below (to avoid pileup).
        Δr = 1e-3
        r_grid = np.arange(0, r_pu+Δr, Δr)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        norm_r = trapz(prob_r, r_grid)
        prob_r /= norm_r
        cdf_r = np.append(0, np.cumsum(prob_r)/np.max(np.cumsum(prob_r)))
        func = interp1d(cdf_r, np.append(0, r_grid))
        r_samples = func(np.random.uniform(size=N_0+N_1+N_2))

        np.testing.assert_almost_equal(trapz(prob_r, r_grid), 1)

        if debugging:
            # The below plot is a sanity check that I ran one time.
            import matplotlib.pyplot as plt
            plt.close('all')
            bin_width = 0.5
            r_bins = np.arange(0,22.5+bin_width,bin_width)
            plt.hist(r_samples, bins=r_bins, normed=False, alpha=0.5,
                    histtype='stepfilled', color='steelblue', edgecolor='none')

            prob_r_on_bins = np.minimum( r_bins**δ, r_pl**δ )
            prob_r_on_bins /= trapz(prob_r_on_bins, r_bins)
            analytic = prob_r_on_bins*(N_0+N_1+N_2)*bin_width
            plt.plot(r_bins, analytic)
            plt.yscale('log')
            plt.savefig('temp3.pdf')

        df['r'] = r_samples
        df.loc[df['has_planet'] == False, 'r'] = np.nan

    elif model_number == 4:
        r_pl = 2 # [Rearth]. Lower and upper bound for truncation.

        # Inverse transform sample to get radii. Drawing from powerlaw
        # distribution above r_pl, and constant below (to avoid pileup).
        Δr = 1e-3
        r_grid = np.arange(0, r_pu+Δr, Δr)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        prob_r[(r_grid > 1.5) & (r_grid < 2)] = 0
        prob_r /= trapz(prob_r, r_grid)
        cdf_r = np.append(0, np.cumsum(prob_r)/np.max(np.cumsum(prob_r)))
        func = interp1d(cdf_r, np.append(0, r_grid))
        r_samples = func(np.random.uniform(size=N_0+N_1+N_2))

        df['r'] = r_samples
        df.loc[df['has_planet'] == False, 'r'] = np.nan

    ###########################################################################
    # Detected planets are those that transit, and whose hosts are searchable #
    ###########################################################################
    Q_g0 = 0.6 # arbitrary geometric transit probability around singles
    Q_g1 = Q_g0
    Q_g2 = Q_g0 * q**(2/3)

    transits_0 = ( np.random.rand((N_0)) < Q_g0 )
    transits_1 = ( np.random.rand((N_1)) < Q_g1 )
    transits_2 = ( np.random.rand((N_2)) < Q_g2 )

    planet_transits = np.concatenate(
        (transits_0, transits_1, transits_2)
        )
    df['planet_transits'] = planet_transits
    # Cleanup: only stars with planets have "planets that transit"
    df.loc[df['has_planet'] == False, 'planet_transits'] = np.nan

    # Assign detected planets
    detected_planets = has_planet.astype(bool) & \
                       planet_transits.astype(bool) & \
                       is_searchable.astype(bool)
    df['has_det_planet'] = detected_planets

    return has_planet, planet_transits, df, Q_g0, r_pu


def observe_planets(df, α):
    '''
    WHAT THE OBSERVER SEES

    Compute apparent radii for detected planets.
    NB. we're saying the detectability is entirely determined by the mass
    ratio (and so the "apparent radius" happens post-detection)
    '''
    df['r_a'] = np.nan

    locind = (df['star_type']=='single') & (df['has_det_planet']==True)
    df.loc[locind, 'r_a'] = df.loc[locind, 'r']

    locind = (df['star_type']=='primary') & (df['has_det_planet']==True)
    this_q = df.loc[locind, 'q']
    this_r = df.loc[locind, 'r']
    df.loc[locind, 'r_a'] = this_r *(1+this_q**α)**(-1/2)

    locind = (df['star_type']=='secondary') & (df['has_det_planet']==True)
    this_q = df.loc[locind, 'q']
    this_r = df.loc[locind, 'r']
    df.loc[locind, 'r_a'] = this_r *(1+this_q**(-α))**(-1/2) * (1/this_q)

    # Assert: only detected planets should have an apparent radius.
    for type_i in ['single', 'primary', 'secondary']:
        assert np.all(np.isfinite(np.array(
            df[(df['star_type']==type_i) & \
            (df['has_det_planet']==True)]['r_a'])))
        assert np.all(~np.isfinite(np.array(
            df[(df['star_type']==type_i) & \
            (df['has_det_planet']==False)]['r_a'])))

    return df


def calculate_true_and_apparent_rates(
        model_number, df, Q_g0, N_0, N_1, N_2, Z_0, Z_1, Z_2, r_pu,
        quickrun, slowrun
        ):
    '''
    Compute the true and apparent rates over bins in true and apparent planet
    radius.

    TRUE RATE DENSITY: defined here to mean true rate density for the stars
    that were _selected_. Note that this is general different from the true
    rate density for a volume-limited sample, or that for singles.

    All of these are the same if all Z_i's are equal.

    APPARENT RATE DENSITY: defined as in 2017/12/05, Kento Masuda's note.

    --------------------

    RETURNS:
        true_dict: true rate of bins in true planet radius
        inferred_dict: apparent rate in bins of apparent planet radius
        outdf: a pandas dataframe with the above, the radius grids, and a few
               other tidbits.

    Also saves `outdf` to the following path:
        '../data/numerics/results_model_{:d}_Zsub2_{:.2f}'.format(
                model_number, Z_2)
    '''

    if model_number == 1 or model_number == 2:
        Δr = 0.01
        radius_bins = np.arange(0, 1+Δr, Δr)
    elif model_number == 3 or model_number == 5:
        Δr = 0.5
        radius_bins = np.arange(0, r_pu+Δr, Δr)
    elif model_number == 4:
        Δr = 0.25
        radius_bins = np.arange(0, r_pu+Δr, Δr)
        #radius_bins = np.logspace(-2,4/3,21) # thought about it, opted against

    true_dict = {}
    inferred_dict = {}
    N_p_at_r_p_inferred = 0
    N_p_at_r_p_true = 0

    types = ['single','primary','secondary']
    for type_i in types:

        # True rate densities
        r = df[(df['star_type']==type_i) & (df['has_planet'])]['r']
        N_p, edges = np.histogram(r, bins=radius_bins)

        true_dict[type_i] = {}
        true_dict[type_i]['N_p'] = N_p
        true_dict[type_i]['edges'] = edges

        # Inferred rate densities
        r_a = df[(df['star_type']==type_i) & (df['has_det_planet'])]['r_a']
        N_p_det, edges = np.histogram( r_a, bins=radius_bins )

        inferred_dict[type_i] = {}
        inferred_dict[type_i]['N_p'] = N_p_det/Q_g0
        inferred_dict[type_i]['edges'] = edges

        if model_number == 1 or model_number == 2:
            rs = np.array(df['r'])
            assert np.all(rs[np.greater(rs, np.zeros_like(rs))] == rs[0])
            r_p = rs[0]
            N_p_at_r_p_inferred += (len(r_a[r_a == r_p])/Q_g0)
            N_p_at_r_p_true += len(r[r == r_p])

    N_tot = N_0 + N_1 + N_2
    true_dict['Λ'] = (true_dict['single']['N_p'] + \
                     true_dict['primary']['N_p'] + \
                     true_dict['secondary']['N_p'])/N_tot
    true_dict['r'] = radius_bins

    N_tot_apparent = N_0 + N_1
    inferred_dict['Λ'] = \
                     (inferred_dict['single']['N_p'] + \
                     inferred_dict['primary']['N_p'] + \
                     inferred_dict['secondary']['N_p'])/N_tot_apparent
    inferred_dict['r'] = radius_bins

    outdf = pd.DataFrame(
            {'bin_left': edges[:-1],
             'true_Λ': true_dict['Λ'],
             'inferred_Λ': inferred_dict['Λ'],
             'true_single_Λ': true_dict['single']['N_p']/N_0,
             'true_primary_Λ': true_dict['primary']['N_p']/N_1,
             'true_secondary_Λ': true_dict['secondary']['N_p']/N_2,
             'frac_inferred_single_Λ': inferred_dict['single']['N_p']/N_tot_apparent,
             'frac_inferred_primary_Λ': inferred_dict['primary']['N_p']/N_tot_apparent,
             'frac_inferred_secondary_Λ': inferred_dict['secondary']['N_p']/N_tot_apparent
            }
            )

    savdir = '../data/numerics/'
    fname = 'results_model_{:d}_Zsub2_{:.2f}_rpu_{:.1f}'.format(
            model_number, Z_2, r_pu)
    if quickrun:
        fname = 'quickrun_' + fname
    outdf.to_csv(savdir+fname+'.out', index=False)
    print('wrote output to {:s}'.format(savdir+fname+'.out'))

    # The sum of the true distribution's histogrammed rate densities is the true
    # average occurrence rate, to three decimals.
    if slowrun:
        np.testing.assert_almost_equal(
                np.sum(outdf['true_Λ']),
                (N_0*Z_0 + N_1*Z_1 + N_2*Z_2)/N_tot,
                decimal=3
                )

    if model_number == 1 or model_number == 2:

        return true_dict, inferred_dict, N_p_at_r_p_true, N_p_at_r_p_inferred

    else:

        return true_dict, inferred_dict, np.nan, np.nan


def run_unit_tests(
        model_number, true_dict, inferred_dict, df,
        N_0, N_1, N_2,
        N_p_at_r_p_true, N_p_at_r_p_inferred, r_pu,
        α, β, γ, δ,
        Z_0, Z_1, Z_2, BF,
        debugging=False
        ):
    '''
    Check numerical results against analytic predictions.
    '''

    N_tot = N_0 + N_1 + N_2
    N_tot_apparent = N_0 + N_1
    μ = N_1/N_0

    if model_number == 1 or model_number== 2:

        Λ_inferred = N_p_at_r_p_inferred / N_tot_apparent
        Λ_true = N_p_at_r_p_true / N_tot
        X_Λ_rp_numerics =  Λ_inferred / Λ_true


    if model_number == 1:

        # TEST: Analytic vs numerical value at the true planet radius, r_p, for
        # model #1. 

        w_a = (1+μ)**(-1)
        w_b = μ*(1+μ)**(-1)

        w0 = (1+2*μ)**(-1)
        w1 = μ * (1+2*μ)**(-1)
        w2 = μ * (1+2*μ)**(-1)

        X_Λ_rp_analytic = w_a/(w0+w1+w2)
        print(X_Λ_rp_analytic)

        # when using larger N_0, gets down to 3 or 4 decimals
        np.testing.assert_almost_equal(
                X_Λ_rp_numerics,
                X_Λ_rp_analytic,
                decimal=2
                )


    elif model_number == 2:

        # TEST #1: Verify that the inferred rate value at the true planet
        # radius agrees with analytic expectation.

        c_0 = 0.494087 # mathematica, 171011_integrals.nb
        X_Λ_rp_analytic = 3*c_0*Z_0/(Z_0+Z_1+Z_2)

        np.testing.assert_almost_equal(
                X_Λ_rp_numerics,
                X_Λ_rp_analytic,
                decimal=2
                )

        # TEST #2: Verify that the inferred rate from 0.73r_true to 0.99r_true
        # is that predicted analytically, to within 0.5%. (No more precise b/c
        # of Poisson noise in bins).

        q_grid = np.arange(1e-6, 1+1e-6, 1e-6)
        prob_ml_q = q_grid**β * (1+q_grid**α)**(3/2)
        norm_q = trapz(prob_ml_q, q_grid)

        r_a_anal = np.array(inferred_dict['r'])[1:] - \
                   np.diff(inferred_dict['r'])/2
        r_anal = np.array(true_dict['r'])[1:] - \
                 np.diff(true_dict['r'])/2
        np.testing.assert_array_equal(r_anal, r_a_anal)

        assert np.allclose(
                np.diff(inferred_dict['r']), np.diff(inferred_dict['r'])[0])

        bin_width = np.diff(inferred_dict['r'])[0]

        rs = np.array(df['r'])
        assert np.all(rs[np.greater(rs, np.zeros_like(rs))] == rs[0])
        r_p = rs[0]

        # Derived 2017/11/30.1. Has the same functional form as the one derived
        # on 2017/10/10.1 (but has explicit normalization).
        I_1 = Z_1/norm_q * 2/α * r_p/(r_a_anal**2) * \
                              ( (r_p/r_a_anal)**2 -1 )**((β-α+1)/α) * \
                              (r_p/r_a_anal)**4
        I_1[r_a_anal == 0] = 0
        I_1[r_a_anal == 1] = 0

        μ = N_1/N_0
        Λ_primaries_analytic = bin_width * μ/(1+μ) * I_1

        # Expect everything from r_p/sqrt(2) to just below r_p to be from the
        # primaries. Going down to 0.73, not 1/sqrt(2) because the funky r_a(q)
        # relation for the i=2 case goes ABOVE sqrt(2) (!). This is lets us
        # verify the numerics work regardless.
        vals_num = inferred_dict['Λ'][(r_a_anal>0.73) & (r_a_anal<0.99)]
        vals_anal = Λ_primaries_analytic[(r_a_anal>0.73) & (r_a_anal<0.99)]

        # Assert the numeric and analytic distributions agree to 0.5%.
        assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=5e-3)

        if debugging:
            # A one-time sanity check.
            plt.close('all')
            f,ax = plt.subplots()
            ax.plot(r_anal[(r_anal>0.73) & (r_anal<0.99)],
                    vals_num,
                    label='numeric')
            ax.plot(r_anal[(r_anal>0.73) & (r_anal<0.99)],
                    vals_anal,
                    label='analytic')
            ax.set_yscale('log')
            ax.legend(loc='best')
            f.savefig('temp0.pdf')

        del vals_num, vals_anal


    elif model_number == 3:

        r_a_anal = np.array(inferred_dict['r'])[1:] - np.diff(inferred_dict['r'])/2
        r_anal = np.array(true_dict['r'])[1:] - np.diff(true_dict['r'])/2

        # Compare the analytic and numeric (true) Λ(r) distributions, above
        # 2r_\oplus
        Λ_num = np.array(true_dict['Λ'])

        r_pl = 2 # [Rearth]. Lower and upper bound for truncation.
        Δr = 1e-3
        r_grid = np.arange(0, r_pu+Δr, Δr)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        norm_r = trapz(prob_r, r_grid)

        Λ_anal = r_anal**δ/norm_r * (
                Z_0/(1+2*μ) + (Z_1+Z_2) * μ/(1+2*μ)
                )

        vals_num = Λ_num[(r_anal>3) & (r_anal<=22)]

        assert np.allclose(np.diff(true_dict['r']), np.diff(true_dict['r'])[0])
        bin_width = np.diff(true_dict['r'])[0]
        vals_anal = bin_width * Λ_anal[(r_anal>3) & (r_anal<=22)]

        # TEST: ensure that the numerically realized values for the true
        # distribution are within 1% of the expected analytic values.
        # (Not exact because of Poisson noise). Print a warning if it seems
        # that Poisson noise is too high.
        assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=1e-2)
        if not np.max(np.abs(1-vals_num/vals_anal))<0.1:
            print('\n'*10)
            print('WARNING: you should use more points to lower Poisson noise')
            print('\n'*10)

        if debugging:
            import matplotlib.pyplot as plt
            plt.close('all')
            f,ax = plt.subplots()
            ax.plot(r_anal[(r_anal>3) & (r_anal<=22)], vals_num, label='numeric')
            ax.plot(r_anal[(r_anal>3) & (r_anal<=22)], vals_anal, label='analytic')
            ax.set_yscale('log')
            ax.legend(loc='best')
            f.savefig('temp_true.pdf')

        # NOTE: It would be ideal if we could compare the analytic and numeric
        # (apparent) Λ_a(r_a) distributions, above r_a ~= 2r_\oplus and below
        # r_pu/sqrt(2).  By analogy with the fixed planet radius, variable mass
        # ratio model -- this might seemingly work.  However if you think about
        # it more closely, it won't. You can derive an analytic result for the
        # f(r)~r^δ case.  But numerically, we have finite edge effects. A
        # r=20rearth planet can be diluted to much below 20rearth/sqrt(2), if
        # it orbits a secondary.  So to make a sufficiently precise analytic
        # prediction, you need to include the finite edge effects.  This might
        # be possible, but it'd certainly be painful.  We've validated enough
        # other cases that the numerics should be trusted (provided they're
        # tested for a variety of upper cutoffs on the radius distribution).


    elif model_number == 5:

        r_a_anal = np.array(inferred_dict['r'])[1:] - \
                   np.diff(inferred_dict['r'])/2
        r_anal = np.array(true_dict['r'])[1:] - \
                 np.diff(true_dict['r'])/2

        # Compare the analytic and numeric (true) Λ(r) distributions, above
        # 2r_\oplus

        Λ_num = np.array(true_dict['Λ'])

        Δr = 1e-3
        r_pl = 2 # [Rearth]. Lower and upper bound for truncation.
        r_grid = np.arange(0, r_pu+Δr, Δr)
        prob_r = np.minimum( r_grid**δ, r_pl**δ )
        norm_r = trapz(prob_r, r_grid)
        prob_r /= norm_r

        Λ_anal = r_anal**δ/norm_r * (
                Z_0/(1+2*μ) + (Z_1+Z_2) * μ/(1+2*μ)
                )

        vals_num = Λ_num[(r_anal>3) & (r_anal<=22)]

        assert np.allclose(np.diff(true_dict['r']), np.diff(true_dict['r'])[0])
        bin_width = np.diff(true_dict['r'])[0]
        vals_anal = bin_width * Λ_anal[(r_anal>3) & (r_anal<=22)]

        if debugging:
            # The following was a sanity check that I ran one time.
            import matplotlib.pyplot as plt
            plt.close('all')
            f,ax = plt.subplots()
            ax.plot(r_anal[(r_anal>3) & (r_anal<=22)], vals_num/vals_anal,
                    label='numeric/analytic')
            ax.legend(loc='best')
            f.savefig('temp2.pdf')

        # TEST: ensure that the numerically realized values for the true
        # distribution are within 1% of the expected analytic values.  (Not
        # exact because of Poisson noise). Print a warning if it seems that
        # Poisson noise is too high.

        assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=1e-2)
        if not np.max(np.abs(1-vals_num/vals_anal))<0.1:
            print('\n'*10)
            print('WARNING: you should use more points to lower Poisson noise')
            print('\n'*10)

        # Compare the analytic and numeric (apparent) Λ_a(r_a) distributions,
        # above r_a ~= 2r_\oplus and below r_pu/sqrt(2).  (I didn't derive the
        # analytic form from r_pu/sqrt(2) to r_pu, but you expect to see a drop
        # in the numeric value because of those HJs being diluted, but it's a
        # small effect and the resulting numerics look fine, CF. the "temp"
        # test plots below).

        del vals_num, vals_anal
        Λ_a_num = np.array(inferred_dict['Λ'])

        # As derived 2017/11/29.3
        Λ_a_anal = r_a_anal**δ / (norm_r *(1+μ)) * (
                Z_0 + 2**((δ+1)/2) * μ * (Z_1+Z_2)
                )

        vals_num = Λ_a_num[(r_a_anal>3) & (r_a_anal<=r_pu/2**(1/2))]

        assert np.allclose(np.diff(inferred_dict['r']),
                           np.diff(inferred_dict['r'])[0])
        bin_width = np.diff(true_dict['r'])[0]

        vals_anal = bin_width * Λ_a_anal[(r_a_anal>3) & (r_a_anal<=r_pu/2**(1/2))]

        # TEST: ensure that the numerically realized values for the apparent
        # distribution are within 1% of the expected analytic values.

        assert np.isclose(np.mean(vals_num/vals_anal), 1, rtol=1e-2)

        if debugging:
            import matplotlib.pyplot as plt
            plt.close('all')
            f,ax = plt.subplots()
            ax.plot(r_anal[(r_anal>3) & (r_anal<=r_pu/2**(1/2))], vals_num/vals_anal,
                    label='numeric/analytic')
            ax.legend(loc='best')
            f.savefig('temp2.pdf')

            plt.close('all')
            f,ax = plt.subplots()
            ax.plot(r_anal[(r_anal>3) & (r_anal<=r_pu/2**(1/2))], vals_num, label='numeric')
            ax.plot(r_anal[(r_anal>3) & (r_anal<=r_pu/2**(1/2))], vals_anal, label='analytic')
            ax.set_yscale('log')
            ax.legend(loc='best')
            f.savefig('temp1.pdf')

        # TEST: Eq 55 of the 2017/12/04 email to JNW and KM. In the regime for
        # which "power law" applies and there are no edge effects, ensure the
        # analytics and numerics agree to 0.5% relative tolerance. NOTE for
        # this to actually be convincing, it is better to use e.g., δ=-2, or
        # δ=-1.  (Because δ=-2.92 leads to a very small shift). The test passes
        # for these cases. E.g., here μ=0.3142696, if δ=-1, Γ_a/Γ_sel is
        # approximately 1.72 (analytically). The numerics agree.
        X_num = Λ_a_num/Λ_num
        mu = N_1/N_0
        assert np.isclose(
                np.mean(X_num[(r_anal>3) & (r_anal<=r_pu/2**(1/2))]),
                (1+2**((δ+3)/2)*mu)/(1+mu),
                rtol=5e-3
                )


def numerical_transit_survey(
        quickrun,
        model_number,
        Z_2,
        BF,
        r_pu,
        α=3.5,
        β=0,
        γ=0,
        δ=-2.92,
        Z_0=0.5,
        Z_1=0.5
        ):
    '''
    Inputs:
        quickrun = False    # True if you want models to run fast.
        model_number = 3

        BF = 0.44   # binary fraction. BF = n_d / (n_s+n_d). Raghavan+ 2010 solar.

        r_pu = 22.5 # maximum allowed radius (only applicable for the power law
                    # radius distribution)

        α = 3.5     # exponent in L~M^α
        β = 0       # exponent in p_vol_limited(q) ~ q^β
        γ = 0       # exponent in p(has planet|secondary w/ q) ~ q^γ
        δ = -2.92   # exponent in p(r) ~ r^δ, from Howard+ 2012

        Z_0 = 0.5 # fraction of selected singles with planet
        Z_1 = 0.5 # fraction of selected primaries with planet
        Z_2 = 0.5 # fraction of selected secondaries with planet.

    Outputs:
        Four data files (in ../data/):

        'results_model_{:d}_Zsub2_{:.2f}'.format( model_number, Z_2)

        'results_model_{:d}_error_case_{:d}_Zsub2_{:.2f}.out'.format(
            model_number, error_case_number, Z_2)

        The output data files contain binned occurrence rates as a function of
        true and apparent radius.

    '''

    print('Running Model #{:d}.'.format(model_number))

    debugging = False

    slowrun = not quickrun

    check_inputs(γ, Z_2, model_number, BF, r_pu)

    N_0, N_1, N_2, is_searchable, df, q = \
            make_stellar_population(quickrun, model_number, BF, α, β, γ, δ)

    has_planet, planet_transits, df, Q_g0, r_pu = \
            make_planet_population(model_number, df, q, is_searchable, Z_0,
                    Z_1, Z_2, N_0, N_1, N_2, δ, r_pu, debugging=debugging)

    df = observe_planets(df, α)

    true_dict, inferred_dict, N_p_at_r_p_true, N_p_at_r_p_inferred = \
            calculate_true_and_apparent_rates(
            model_number, df, Q_g0, N_0, N_1, N_2, Z_0, Z_1, Z_2, r_pu,
            quickrun, slowrun
            )

    run_unit_tests(
            model_number, true_dict, inferred_dict, df,
            N_0, N_1, N_2,
            N_p_at_r_p_true, N_p_at_r_p_inferred, r_pu,
            α, β, γ, δ,
            Z_0, Z_1, Z_2, BF, debugging=debugging
            )


if __name__ == '__main__':

    np.random.seed(42)

    parser = argparse.ArgumentParser(
        description='Run a numerical transit survey.')

    parser.add_argument('-qr', '--quickrun', action='store_true',
        help='use --quickrun if you want models to run fast.')
    parser.add_argument('-mn', '--modelnumber', type=int, default=None,
        help='1, 2 or 3.')
    parser.add_argument('-ztwo', '--ZsubTwo', type=float, default=None,
        help='integrated occ rate for secondaries')
    parser.add_argument('-bf', '--binaryfrac', type=float, default=None,
        help='BF = n_b/(n_s+n_b), for n number density')
    parser.add_argument('-rpu', '--upperradiuscutoff', type=float, default=None,
        help='the maximum allowed planet radius in units of Rearth')

    args = parser.parse_args()

    numerical_transit_survey(
        args.quickrun,
        args.modelnumber,
        args.ZsubTwo,
        args.binaryfrac,
        args.upperradiuscutoff)
