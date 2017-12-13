# -*- coding: utf-8 -*-
from __future__ import division, print_function
import matplotlib as mpl
mpl.use('pgf')
pgf_with_custom_preamble = {
    'pgf.texsystem': 'pdflatex', # xelatex is default; i don't have it
    'font.family': 'serif', # use serif/main font for text elements
    'text.usetex': True,    # use inline math for ticks
    'pgf.rcfonts': False,   # don't setup fonts from rc parameters
    'pgf.preamble': [
        '\\usepackage{amsmath}',
        '\\usepackage{amssymb}'
        ]
    }
mpl.rcParams.update(pgf_with_custom_preamble)

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

    plt.close('all')

    # Make summary plot
    fname = '../data/numerics/integrate_model_{:d}_Zsub2_{:.2f}.out'.format(
            model_number, Z_2)

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

    if model_number in [1,2]:
        xvals = np.array(df['r'])
        yinferred = np.array(df['Γ_a'])
        xsingle = np.append(xvals,1)
        ytruesingle = np.append(np.array(df['Γ_0']), 0)
    elif model_number in [3,4,5,6]:
        xvals = np.array(df['r'])
        yinferred = np.array(df['Γ_a'])
        xsingle = xvals
        ytruesingle = np.array(df['Γ_0'])

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    ax.plot(xvals, yinferred, label='inferred', c=colors[0])

    if standardlines:

        ax.plot(xsingle, ytruesingle, label='true (single)', linestyle='-',
                c=colors[1])

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
        ax.legend(loc='upper left',fontsize='medium')
    else:
        ax.legend(loc='best',fontsize='medium')

    ax.set_xlabel('planet radius, $r$ [$r_\oplus$]', fontsize='large')

    ax.set_ylabel('rate density, $\Gamma(r)$ [$r_\oplus^{-1}$]', fontsize='large')

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    if logy and not logx:
        outname = '../results/int_rate_density_vs_radius_logy_model_'+\
                 repr(model_number)
    elif not logy and logx:
        outname = '../results/int_rate_density_vs_radius_logx_model_'+\
                 repr(model_number)
    elif logx and logy:
        outname = '../results/int_rate_density_vs_radius_logs_model_'+\
                 repr(model_number)
    else:
        outname = '../results/int_rate_density_vs_radius_model_'+\
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
            inds = df['r'] > lower_bound
            Λ_HJ_true_sing = trapz(df[inds]['Γ_0'], df[inds]['r'])
            Λ_HJ_inferred = trapz(df[inds]['Γ_a'], df[inds]['r'])
            #Λ_HJ_true_sel = np.sum(df[inds]['true_Λ'])
            #Λ_HJ_true_sing = np.sum(df[inds]['true_single_Λ'])
            #Λ_HJ_inferred = np.sum(df[inds]['inferred_Λ'])

            Λ_true_sing = trapz(df['Γ_0'], df['r'])
            Λ_inferred = trapz(df['Γ_a'], df['r'])
            #Λ_true_sel = np.sum(df['true_Λ'])
            #Λ_true_sing = np.sum(df['true_single_Λ'])
            #Λ_inferred = np.sum(df['inferred_Λ'])

            txt = \
            '''
            with $r>${:.1f}$R_\oplus$,
            single $\Lambda$_HJ_true:      {:.4f} planets per star
            $\Lambda$_HJ_inferred:  {:.4f} planets per star
            true(single)/inferred:  {:.2f}.

            Integrated over all $r$,
            single $\Lambda$_true:           {:.3f} planets per star
            $\Lambda$_inferred:     {:.3f} planets per star
            true(single)/inferred:  {:.2f}.
            '''.format(
            lower_bound,
            Λ_HJ_true_sing,
            Λ_HJ_inferred,
            Λ_HJ_true_sing/Λ_HJ_inferred,
            Λ_true_sing,
            Λ_inferred,
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
            Z_0 = 0.5
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

    make_plot(2, logy=True, Z_2=0.5)
    make_plot(5, logx=True, logy=True, Z_2=0.5)
    make_plot(6, logx=True, logy=True, Z_2=0.5)

    make_plot(3, Z_2=0.5)
    make_plot(3, withtext=True, Z_2=0.5)
    make_plot(3, logx=False, logy=True, Z_2=0.5)
    make_plot(3, logx=True, logy=True, Z_2=0.5)

    make_plot(4, Z_2=0.5, xcut=True)
    make_plot(4, logy=True, Z_2=0.5, xcut=True)
