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
        mainly for filenaming
    '''

    assert Z_2 > -1
    Z_0 = 0.5

    plt.close('all')

    if not brokenx:
        f, ax = plt.subplots(figsize=(4,4))
    elif many_Zs:
        f, ax = plt.subplots(figsize=(5,4))
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
        raise NotImplementedError
    elif model_number in [3,4,5,6,7]:
        xvals = np.array(df['r'])
        y_0= np.array(df['ndet_a0'])
        y_1= np.array(df['ndet_a1'])
        y_2= np.array(df['ndet_a2'])
        y_a = y_0 + y_1 + y_2

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    c = ['#4ca8e8', '#85c5f2'] # shades of blue for vs function of Z_2

    ndetstr = '/$n_{\mathrm{det}}^{\mathrm{tot}}(r_{\mathrm{a}})$'
    if not many_Zs:
        #ax.plot(xvals, y_a, label=ndetstr, c=colors[0])
        ax.plot(xvals, y_0/y_a,
                label='$n_{\mathrm{det}}^0(r_{\mathrm{a}})$'+ndetstr,
                c=colors[0], linestyle='-')
        ax.plot(xvals, y_1/y_a,
                label='$n_{\mathrm{det}}^1(r_{\mathrm{a}})$'+ndetstr,
                c=colors[0], linestyle='--')
        ax.plot(xvals, y_2/y_a,
                label='$n_{\mathrm{det}}^2(r_{\mathrm{a}})$'+ndetstr,
                c=colors[0], linestyle=':')
    else:
        for ind, Z_2 in enumerate([0.5,0.25,0]):
            fname = '../data/numerics/integrate_model_{:d}_Zsub2_{:.2f}_rpu_{:.2f}.out'.format(
                    model_number, Z_2, r_pu)

            df = pd.read_csv(fname)

            xvals = np.array(df['r'])
            y_0= np.array(df['ndet_a0'])
            y_1= np.array(df['ndet_a1'])
            y_2= np.array(df['ndet_a2'])
            y_a = y_0 + y_1 + y_2

            c = ['#4ca8e8', '#85c5f2'] # shades of blue for vs function of Z_2
            ls = ['--', ':'] # shades of blue for vs function of Z_2

            Zstr = ', $Z_2/Z_0$={:.1f}'.format(Z_2/Z_0)
            ax.plot(xvals, y_0/y_a,
                    label='$n_{\mathrm{det}}^0(r_{\mathrm{a}})$'+ndetstr+Zstr,
                    c=colors[ind], linestyle='-')
            ax.plot(xvals, y_1/y_a,
                    label='$n_{\mathrm{det}}^1(r_{\mathrm{a}})$'+ndetstr+Zstr,
                    c=colors[ind], linestyle='--')
            ax.plot(xvals, y_2/y_a,
                    label='$n_{\mathrm{det}}^2(r_{\mathrm{a}})$'+ndetstr+Zstr,
                    c=colors[ind], linestyle=':')

            print(Z_2, min(y_1/y_a), max(y_1/y_a), max(y_1/y_a)/min(y_1/y_a))


    if many_Zs:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 5/5, box.height])
	# put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if model_number != 7 and not many_Zs:
        ax.legend(loc='best',fontsize='medium')

    if model_number ==7 and not many_Zs:
        ax.legend(loc='upper left',fontsize='small')

    ax.set_xlabel('apparent planet radius, $r_{\mathrm{a}}$ [$r_\oplus$]', fontsize='large')

    #ax.set_ylabel('number of planets detected [$N_{\mathrm{s}}^0(r_{\mathrm{a}})\cdot p_{\mathrm{tra}}$]', fontsize='large')
    ax.set_ylabel('fraction of detected planets', fontsize='large')

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    if logy and not logx:
        outname = '../results/ndet_vs_radius_logy_model_'+\
                 repr(model_number)
    elif not logy and logx:
        outname = '../results/ndet_vs_radius_logx_model_'+\
                 repr(model_number)
    elif logx and logy:
        outname = '../results/ndet_vs_radius_logs_model_'+\
                 repr(model_number)
    else:
        outname = '../results/ndet_vs_radius_model_'+\
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

    if many_Zs and logx:
        ax.set_xlim([0.2,22])

    if model_number == 3 or model_number==7:

        #Z_2 = 0.5
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

    #import IPython; IPython.embed()
    f.savefig(outname+'.pdf', bbox_inches='tight')
    print('saved {:s}'.format(outname+'.pdf'))


if __name__ == '__main__':

    #note: not running
    #make_plot(2, logy=True, Z_2=0.5, r_pu=1)

    #make_plot(5, logx=True, logy=True, Z_2=0.5, r_pu=22.5, standardlines=False)
    #make_plot(6, logx=True, logy=True, Z_2=0.5, r_pu=22.5, standardlines=False)

    #make_plot(3, Z_2=0.5, r_pu=22.5, standardlines=False)
    #make_plot(3, withtext=True, Z_2=0.5, r_pu=22.5, standardlines=False)
    #make_plot(3, logx=False, logy=True, Z_2=0.5, r_pu=22.5, standardlines=False)
    #make_plot(3, logx=True, logy=True, Z_2=0.5, r_pu=22.5, standardlines=False)

    #make_plot(4, Z_2=0.5, xcut=True, r_pu=22.5, standardlines=False)
    #make_plot(4, logy=True, Z_2=0.5, xcut=True, r_pu=22.5, standardlines=False)

    ## Change as a function of Z_2/Z_0
    #for Z_2 in [0, 0.25]:
    #    make_plot(3, Z_2=Z_2, r_pu=22.5, standardlines=False)
    #    make_plot(3, logy=True, Z_2=Z_2, r_pu=22.5, standardlines=False)
    #    make_plot(3, logx=True, logy=True, Z_2=Z_2, r_pu=22.5, standardlines=False)
    #    make_plot(3, withtext=True, stdout=True, Z_2=Z_2, r_pu=22.5, standardlines=False)

    # As a function of Z_2/Z_0, on the same plot (...)
    make_plot(3, r_pu=22.5, many_Zs=True, Z_2=0.5, standardlines=False,
            logx=True)
    #make_plot(3, r_pu=22.5, many_Zs=True, Z_2=0.5, standardlines=False)
    #make_plot(4, r_pu=22.5, many_Zs=True, Z_2=0.5, xcut=True, standardlines=False)
    make_plot(4, r_pu=22.5, many_Zs=True, Z_2=0.5, xcut=True,
            standardlines=False, logx=True)

    make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5, xcut=True,
            standardlines=False, logx=True)
    #make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5, standardlines=False)
    #make_plot(7, r_pu=22.5, many_Zs=True, Z_2=0.5, withtext=True, standardlines=False)

