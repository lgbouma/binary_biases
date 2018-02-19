# -*- coding: utf-8 -*-
'''
plot the mass ratio distribution for a volume-limited sample, and a
magnitude-limited sample.
'''
from __future__ import division, print_function
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from scipy.integrate import trapz

f,ax = plt.subplots(figsize=(4,4))

q = np.arange(-0.1,1.1+1e-3,1e-3)
β = 0
α = 3.5

yvals = q**β * (1+q**α)**(3/2)
yvals[(q>1)|(q<0)] = 0
ax.plot(q, yvals/trapz(yvals, q), label='magnitude-limited', zorder=2)

yvals = np.ones_like(q)
yvals[(q>1)|(q<0)] = 0
ax.plot(q, yvals, label='volume-limited', zorder=1)

ax.legend(loc='best',fontsize='medium')

ax.set_xlabel('binary mass ratio, $q=M_2/M_1$', fontsize='large')
#ax.set_ylabel('$p(\mathrm{draw\ }q\,|\,\mathrm{system\ is\ binary})$', fontsize='large')
ax.set_ylabel('relative probability', fontsize='large')
ax.set_ylim([-0.1,2.3])
ax.set_xlim([-0.05,1.05])

outname = '../results/mass_ratio_distribution'
f.savefig(outname+'.pdf', bbox_inches='tight')
