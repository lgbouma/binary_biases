'''
visualize the selected and searchable volumes for single stars and twin
binaries.
'''
# -*- coding: utf-8 -*-
from __future__ import division, print_function
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import proj3d

def make_sphere(radius, txt, ax, ind):
    # radius: float
    # txt: what you want annotation to say
    # ax: Axes3D instance
    # ind: deals w/ obnoxious text placement

    # Make data
    u = np.linspace(-np.pi, np.pi, 500)
    v = np.linspace(0, np.pi, 500)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    surf = ax.plot_surface(x, y, z, color='skyblue', alpha=0.3, linewidth=0,
            zorder=1)
    # Make sure it's rasterized to prevent bloated figures. NB. matplotlib
    # throws an error saying that this is ignored, but it's clearly not, b/c it
    # slims the output plot by a factor of 10.
    surf.set_rasterized(True)

    ax.plot(radius*np.cos(u), radius*np.abs(np.sin(u)), 0, color='k', alpha=0.9,
            zorder=2)

    _x, _y = np.linspace(0,radius,100), np.linspace(0,radius,100)
    ax.plot(-np.sqrt(_x**2+_y**2)/2, np.sqrt(_x**2+_y**2)/2, np.zeros_like(_x),
            color='k', linestyle=':', zorder=3)

    minbound, maxbound = -1, 1
    ax.auto_scale_xyz(
            [minbound, maxbound], [minbound, maxbound], [minbound, maxbound])

    ax.view_init(30, 90)

    #######################
    # MAKE THE ANNOTATION #
    #######################
    #end of arrow
    x1, y1, _ = proj3d.proj_transform(-2**(-3/2)*radius, 2**(-3/2)*radius, 0,
                                      ax.get_proj())
    #beginning of arrow
    prefactor = 1.2 if ind != 2 else 1
    x2, y2, _ = proj3d.proj_transform(-prefactor*radius, 1*radius, 1*radius, ax.get_proj())

    ann = ax.annotate(
        txt,
        xy=(x1, y1), xycoords='data',
        xytext=(x2, y2), textcoords='data',
        ha='left', va='bottom',
        fontsize='large',
        bbox = dict(boxstyle='round,pad=0', fc='white', alpha=0, linewidth=0),
        arrowprops = dict(
            arrowstyle='->',
            color='black',
            connectionstyle='angle3,angleA=0,angleB=90'),
        zorder=4)

    ax.set_axis_off()
    ax.patch.set_alpha(0.) # make background transparent


fig = plt.figure(figsize=(4,4))

inds = [1,2,3,4]
radii = [1, 2**(1/2), 1, 2**(-1/2)]
txts = ['$1$', '$\sqrt{2}$', '$1$', '$1/\sqrt{2}$']

for ind, r, txt in zip(inds, radii, txts):

    ax = fig.add_subplot(2,2,ind,projection='3d')

    make_sphere(r, txt, ax, ind)

fig.tight_layout(h_pad=-4, w_pad=-4) #height & width pad

#Overplot the titles in figure coordinate system
ax.text2D(.3, .05, 'single stars', transform=fig.transFigure,
     fontsize='large', va='bottom', ha='center')
ax.text2D(.7, .05, 'twin binaries', transform=fig.transFigure,
     fontsize='large', va='bottom', ha='center')
ax.text2D(.05, .3, 'searchable', transform=fig.transFigure,
     fontsize='large', va='center', ha='left', rotation=90)
ax.text2D(.05, .7, 'selected', transform=fig.transFigure,
        fontsize='large', va='center', ha='left', rotation=90)

fig.savefig('../results/visualize_volumes.pdf', dpi=300,
        bbox_inches='tight', transparent=True)
