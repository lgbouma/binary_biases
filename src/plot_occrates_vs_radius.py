# -*- coding: utf-8 -*-
'''
This script plots the output of the Monte Carlo simulation.

first, for models 1,2,3,4:
>>> python numerical_models.py
then:
>>> python plot_occrates_vs_radius.py

makes plot for Model # 1 of paper
'''
from __future__ import division, print_function
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from brokenaxes import brokenaxes

def make_plot(model_number, logx=None, logy=None, withtext=None,
        stdout=False, brokenx=None, Z_2=None, xcut=None, r_pu=None,
        standardlines=True):
    '''
    withtext: includes summary text from integration on plot

    stdout: write "withtext" text to stdout

    xcut: mainly for model #1, overplotted over different axes

    standardlines:
        if True: plots "inferred", "true (selected)", and "true (single)" rates.
        if False, plots "inferred", "frac (single)", "frac (primary)", "frac
            (secondary)" and "true (single)".
        by default, True
    '''

    assert Z_2 > -1
    Z_0 = 0.5

    plt.close('all')

    # Make summary plot
    fname = '../data/numerics/results_model_{:d}_Zsub2_{:.2f}_rpu_{:.1f}.out'.format(
            model_number, Z_2, r_pu)

    df = pd.read_csv(fname)

    if not brokenx:
        f, ax = plt.subplots(figsize=(4,4))
    else:
        #cf https://github.com/bendichter/brokenaxes/blob/master/brokenaxes.py
        f = plt.figure(figsize=(4,4))
        bigax = brokenaxes(
                xlims=((0.695, .715), (.985, 1.005)),
                d=0.02,
                tilt=85,
                hspace=.05,
                despine=True)
        ax=bigax

    if model_number == 3 or model_number == 4:
        xvals = np.append(0, df['bin_left'])
        ytrueselected = np.append(0, df['true_Λ'])
        ytruesingle = np.append(0, df['true_single_Λ'])
        yinferred = np.append(0, df['inferred_Λ'])
    elif model_number == 1 or model_number == 2:
        xvals = np.append(df['bin_left'],[1,1.1])
        ytrueselected = np.append(df['true_Λ'],[0,0])
        ytruesingle = np.append(df['true_single_Λ'],[0,0])
        yinferred = np.append(df['inferred_Λ'],[0,0])

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    lambdastr = r'$\Lambda_{\mathrm{a}}(r_{\mathrm{a}})$'
    ax.step(xvals, yinferred, where='post', label=lambdastr, c=colors[0])

    if standardlines:

        #ax.step(xvals, ytrueselected, where='post', label='true (selected)',
        #        c=colors[1])

        ax.step(xvals, ytruesingle, where='post', label='$\Lambda(r)$',
                linestyle='-', c=colors[1])

    elif not standardlines:

        yfracsingle = np.append(0, df['frac_inferred_single_Λ'])
        yfracprimary = np.append(0, df['frac_inferred_primary_Λ'])
        yfracsecondary = np.append(0, df['frac_inferred_secondary_Λ'])

        ax.step(xvals, yfracsingle, where='post', label='frac (single)',
                linestyle=':', c=colors[0])
        ax.step(xvals, yfracprimary, where='post', label='frac (primary)',
                linestyle='-.', c=colors[0])
        ax.step(xvals, yfracsecondary, where='post', label='frac (secondary)',
                linestyle='--', c=colors[0])
        ax.step(xvals, ytruesingle, where='post', label='true (single)',
                linestyle='-', c=colors[1])

    if brokenx:
        leg = ax.legend(loc='upper left',fontsize='medium')
        leg.get_frame().set_linewidth(0.)
        leg.get_frame().set_facecolor('none')
    else:
        ax.legend(loc='best',fontsize='medium')

    ax.set_xlabel('planet radius [$R_{\oplus}$]', fontsize='large')

    ax.set_ylabel('number of planets per star', fontsize='large')

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    if logx or logy:
        outname = '../results/occ_rate_vs_radius_logs_model_'+\
                 repr(model_number)
    else:
        outname = '../results/occ_rate_vs_radius_model_'+\
                 repr(model_number)
    if brokenx:
        outname += '_brokenx'

    if not standardlines:
        outname += '_fraclines'

    if not brokenx:
        if (model_number == 1 or model_number == 2) and not (logx or logy):
            ax.set_xlim([0.5,1.02])
        elif (model_number == 1 or model_number == 2) and (logx or logy):
            ax.set_xlim([-0.02,1.02])

    if xcut:
        ax.set_xlim([-0.3,5.3])

    if model_number == 3:
        # Assess HJ rate difference.
        from scipy.integrate import trapz

        #Howard 2012 boundary #1 and boundary #2:
        for lower_bound in [5.6,8]:
            inds = df['bin_left'] > lower_bound
            #Λ_HJ_true = trapz(df[inds]['true_Λ'], df[inds]['bin_left'])
            #Λ_HJ_inferred = trapz(df[inds]['inferred_Λ'], df[inds]['bin_left'])
            Λ_HJ_true_sel = np.sum(df[inds]['true_Λ'])
            Λ_HJ_true_sing = np.sum(df[inds]['true_single_Λ'])
            Λ_HJ_inferred = np.sum(df[inds]['inferred_Λ'])

            #Λ_true = trapz(df['true_Λ'], df['bin_left'])
            #Λ_inferred = trapz(df['inferred_Λ'], df['bin_left'])
            Λ_true_sel = np.sum(df['true_Λ'])
            Λ_true_sing = np.sum(df['true_single_Λ'])
            Λ_inferred = np.sum(df['inferred_Λ'])

            txt = \
            '''
            with $r>${:.1f}$R_\oplus$,
            selected $\Lambda$_HJ_true:      {:.4f} planets per star
            single   $\Lambda$_HJ_true:      {:.4f} planets per star
            $\Lambda$_HJ_inferred:  {:.4f} planets per star
            true(selected)/inferred:  {:.2f}.
            true(single)/inferred:  {:.2f}.

            Integrated over all $r$,
            selected $\Lambda$_true:         {:.3f} planets per star
            single $\Lambda$_true:           {:.3f} planets per star
            $\Lambda$_inferred:     {:.3f} planets per star
            true(selected)/inferred:  {:.2f}.
            true(single)/inferred:  {:.2f}.
            '''.format(
            lower_bound,
            Λ_HJ_true_sel,
            Λ_HJ_true_sing,
            Λ_HJ_inferred,
            Λ_HJ_true_sel/Λ_HJ_inferred,
            Λ_HJ_true_sing/Λ_HJ_inferred,
            Λ_true_sel,
            Λ_true_sing,
            Λ_inferred,
            Λ_true_sel/Λ_inferred,
            Λ_true_sing/Λ_inferred,
            )
            if stdout:
                print(txt)

        if withtext:
            ax.text(0.96,0.5,txt,horizontalalignment='right',
                    verticalalignment='center',
                    transform=ax.transAxes, fontsize='x-small')
            outname += '_withtext'

        else:
            txt = '$Z_2/Z_0 =\ ${:.1f}'.format(Z_2/Z_0)
            ax.text(0.96,0.5,txt,horizontalalignment='right',
                    verticalalignment='center',
                    transform=ax.transAxes, fontsize='x-small')

        if isinstance(Z_2,float) or isinstance(Z_2,int):
            outname += '_Zsub2_{:.2f}'.format(Z_2)

        if isinstance(r_pu,float) or isinstance(r_pu,int):
            outname += '_rpu_{:.1f}'.format(r_pu)

    if xcut:
        outname += '_xcut'

    f.savefig(outname+'.pdf', bbox_inches='tight')



if __name__ == '__main__':

    # for r_pu, model_number in zip([1, 1, 22.5], [1,2,3]):

    #     make_plot(model_number, Z_2=0.5, r_pu=r_pu)

    make_plot(1, brokenx=True, Z_2=0.5, r_pu=1)
    # make_plot(2, logy=True, Z_2=0.5, r_pu=1)
    # make_plot(3, logy=True, Z_2=0.5, r_pu=22.5)
    # make_plot(3, withtext=True, stdout=True, Z_2=0.5, r_pu=22.5)

    # make_plot(4, Z_2=0.5, r_pu=22.5)
    # make_plot(4, xcut=True, Z_2=0.5, r_pu=22.5)

    # # Change as a function of Z_2/Z_0
    # for Z_2 in [0, 0.25]:
    #     make_plot(3, Z_2=Z_2, r_pu=22.5)
    #     make_plot(3, logy=True, Z_2=Z_2, r_pu=22.5)
    #     make_plot(3, withtext=True, stdout=True, Z_2=Z_2, r_pu=22.5)

    # for Z_2 in [0, 0.25, 0.5]:
    #     make_plot(3, Z_2=Z_2, r_pu=22.5, standardlines=False)

    # # Change as a function of r_pu
    # for r_pu in [15,17.5,20,22.5,25]:
    #     make_plot(3, Z_2=0.5, r_pu=r_pu)
    #     make_plot(3, logy=True, Z_2=0.5, r_pu=r_pu)
    #     make_plot(3, withtext=True, stdout=True, Z_2=0.5, r_pu=r_pu)

    # # If you fine-tune both r_pu AND Z_2/Z_0 preferentially, how big of a "HJ
    # # discrepancy" do you get?
    # make_plot(3, withtext=True, stdout=True, Z_2=0., r_pu=15.0)
