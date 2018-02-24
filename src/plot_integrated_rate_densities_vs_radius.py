# -*- coding: utf-8 -*-
from __future__ import division, print_function
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import matplotlib.pyplot as plt, pandas as pd, numpy as np
from scipy.integrate import trapz

from brokenaxes import brokenaxes

def make_plot(model_number, logx=None, logy=None, withtext=None,
        stdout=False, brokenx=None, Z_2=None, xcut=None, r_pu=None,
        standardlines=True, many_Zs=False):
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

    # Make summary plot
    fname = '../data/numerics/integrate_model_{:d}_Zsub2_{:.2f}_rpu_{:.2f}.out'.format(
            model_number, Z_2, r_pu)

    df = pd.read_csv(fname)

    if model_number in [1,2]:
        xvals = np.array(df['r'])
        yinferred = np.array(df['Γ_a'])
        xsingle = np.append(xvals,1)
        ytruesingle = np.append(np.array(df['Γ_0']), 0)
    elif model_number in [3,4,5,6,7]:
        xvals = np.array(df['r'])
        yinferred = np.array(df['Γ_a'])
        xsingle = xvals
        ytruesingle = np.array(df['Γ_0'])

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    c = ['#4ca8e8', '#85c5f2'] # shades of blue for vs function of Z_2

    gammastr = r'$\Gamma_{\mathrm{a}}(r_{\mathrm{a}})$'
    if not many_Zs:
        ax.plot(xvals, yinferred, label=gammastr, c=colors[0])
    else:
        if Z_2/Z_0 == 0.5:
            labelstr = gammastr+', $Z_2/Z_0$={:.1f}'.format(Z_2/Z_0)
        else:
            labelstr = gammastr+', $Z_2/Z_0$={:d}'.format(int(Z_2/Z_0))
        ax.plot(xvals, yinferred, label=labelstr, c=colors[0])

    if many_Zs:
        for ind, Z_2 in enumerate([0.25,0]):
            fname = '../data/numerics/integrate_model_{:d}_Zsub2_{:.2f}_rpu_{:.2f}.out'.format(
                    model_number, Z_2, r_pu)

            df = pd.read_csv(fname)

            xvals = np.array(df['r'])
            yinferred = np.array(df['Γ_a'])

            c = ['#4ca8e8', '#85c5f2'] # shades of blue for vs function of Z_2
            ls = ['--', ':'] # shades of blue for vs function of Z_2

            if Z_2/Z_0 == 0.5:
                labelstr = gammastr+', $Z_2/Z_0$={:.1f}'.format(Z_2/Z_0)
            else:
                labelstr = gammastr+', $Z_2/Z_0$={:d}'.format(int(Z_2/Z_0))

            ax.plot(xvals, yinferred, label=labelstr, c=c[ind],
                    linestyle=ls[ind])


    if standardlines:

        ax.plot(xsingle, ytruesingle, label='$\Gamma_0(r)$', linestyle='-',
                c=colors[1])

    if model_number != 7:
        leg = ax.legend(loc='best',fontsize='medium')
        leg.get_frame().set_linewidth(0.)
        leg.get_frame().set_facecolor('none')
    elif model_number == 7 and not logy:
        leg = ax.legend(loc='upper left',fontsize='small')
        leg.get_frame().set_linewidth(0.)
        leg.get_frame().set_facecolor('none')
    else:
        ax.legend(loc='lower right',fontsize='small')

    ax.set_xlabel('planet radius [$R_\oplus$]', fontsize='large')

    ax.set_ylabel('occurrence rate density [$R_\oplus^{-1}$]', fontsize='large')

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

    if model_number == 3 or model_number==7:

        Z_2 = 0.5
        fname = '../data/numerics/integrate_model_{:d}_Zsub2_{:.2f}_rpu_{:.2f}.out'.format(
                model_number, Z_2, r_pu)
        df = pd.read_csv(fname)

        # Assess HJ rate difference.
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

        elif model_number==3 and not many_Zs:
            txt = '$Z_2/Z_0 =\ ${:.1f}'.format(Z_2/Z_0)
            ax.text(0.96,0.5,txt,horizontalalignment='right',
                    verticalalignment='center',
                    transform=ax.transAxes, fontsize='x-small')

        if isinstance(r_pu,float) or isinstance(r_pu,int):
            outname += '_rpu_{:.1f}'.format(r_pu)

    if (isinstance(Z_2,float) or isinstance(Z_2,int)) and not many_Zs:
        outname += '_Zsub2_{:.2f}'.format(Z_2)

    if xcut:
        outname += '_xcut'

    if many_Zs:
        outname += '_manyZs'

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    f.savefig(outname+'.pdf', bbox_inches='tight')


if __name__ == '__main__':

    #note: not running
    #make_plot(2, logy=True, Z_2=0.5, r_pu=1)

    #make_plot(5, logx=True, logy=True, Z_2=0.5, r_pu=22.5)
    #make_plot(6, logx=True, logy=True, Z_2=0.5, r_pu=22.5)

    #make_plot(3, Z_2=0.5, r_pu=22.5)
    #make_plot(3, withtext=True, Z_2=0.5, r_pu=22.5)
    #make_plot(3, logx=False, logy=True, Z_2=0.5, r_pu=22.5)
    #make_plot(3, logx=True, logy=True, Z_2=0.5, r_pu=22.5)

    #make_plot(4, Z_2=0.5, xcut=True, r_pu=22.5)
    #make_plot(4, logy=True, Z_2=0.5, xcut=True, r_pu=22.5)

    ## Change as a function of Z_2/Z_0
    #for Z_2 in [0, 0.25]:
    #    make_plot(3, Z_2=Z_2, r_pu=22.5)
    #    make_plot(3, logy=True, Z_2=Z_2, r_pu=22.5)
    #    make_plot(3, logx=True, logy=True, Z_2=Z_2, r_pu=22.5)
    #    make_plot(3, withtext=True, stdout=True, Z_2=Z_2, r_pu=22.5)

    ## As a function of Z_2/Z_0, on the same plot (...)
    make_plot(3, r_pu=22.5, many_Zs=True, Z_2=0.5)
    make_plot(4, r_pu=22.5, many_Zs=True, Z_2=0.5, xcut=True)
    make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5)
    #make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5, withtext=True)

    #make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5, logy=True)

    # Change as a function of r_pu
    #for r_pu in [15,17.5,20,22.5,25]:
    #    make_plot(3, Z_2=0.5, r_pu=r_pu)
    #    make_plot(3, logy=True, Z_2=0.5, r_pu=r_pu)
    #    make_plot(3, withtext=True, stdout=True, Z_2=0.5, r_pu=r_pu)

    ## If you fine-tune both r_pu AND Z_2/Z_0 preferentially, how big of a "HJ
    ## discrepancy" do you get?
    #make_plot(3, withtext=True, stdout=True, Z_2=0, r_pu=15)
