import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def main():

    ############################################################
    # Loading data

    Nx = 161
    Ny = 398
    h = 12.5

    l_lim = 0 / 1000
    r_lim = (Ny - 1) * h / 1000
    t_lim = 0.
    b_lim = 2.

    data1 = sio.loadmat('case01.mat')
    data2 = sio.loadmat('case02.mat')
    data3 = sio.loadmat('case03.mat')

    d1 = data1['d'].transpose()
    d2 = data2['d'].transpose()
    d3 = data3['d'].transpose()

    ############################################################
    # Plotting data

    cmap = 'seismic'
    p_lim = 2
    fs1 = 14
    fs2 = 10

    fig = plt.figure(1, figsize=(7, 6))

    ax1 = fig.add_axes((.1, .51, .39, .39))
    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_ticks_position('left')
    ax1.set_yticks(np.linspace(0, 2, 5))
    ax1.yaxis.set_label_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Time (s)')
    ax1.imshow(d1, cmap=cmap, aspect='auto',
               vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.1, 0.05, '(a)', fontsize=fs1, verticalalignment='top')
    plt.text(4.9, 1.9, 'Acoustic', fontsize=fs2, horizontalalignment='right')


    ax2 = fig.add_axes((.51, .51, .39, .39))
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.yaxis.set_ticks_position('right')
    ax2.set_yticks(np.linspace(0, 2, 5))
    ax2.yaxis.set_label_position('right')
    plt.xlabel('Distance (km)')
    plt.ylabel('Time (s)')
    ax2.imshow(d2, cmap=cmap, aspect='auto',
               vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.1, 0.05, '(b)', fontsize=fs1, verticalalignment='top')
    plt.text(4.9, 1.9, 'Viscoacoustic (New)', fontsize=fs2, horizontalalignment='right')

    ax3 = fig.add_axes((.1, .09, .39, .39))
    ax3.xaxis.set_ticks_position('bottom')
    ax3.xaxis.set_label_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.set_yticks(np.linspace(0, 2, 5))
    ax3.yaxis.set_label_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Time (s)')
    ax3.imshow(d3, cmap=cmap, aspect='auto',
               vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.1, 0.05, '(c)', fontsize=fs1, verticalalignment='top')
    plt.text(4.9, 1.9, 'Viscoacoustic (Reference)', fontsize=fs2, horizontalalignment='right')

    ax4 = fig.add_axes((.51, .09, .39, .39))
    ax4.xaxis.set_ticks_position('bottom')
    ax4.xaxis.set_label_position('bottom')
    ax4.yaxis.set_ticks_position('right')
    ax4.set_yticks(np.linspace(0, 2, 5))
    ax4.yaxis.set_label_position('right')
    plt.xlabel('Distance (km)')
    plt.ylabel('Time (s)')
    ax4.imshow(d2-d3, cmap=cmap, aspect='auto',
               vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.1, 0.05, '(d)', fontsize=fs1, verticalalignment='top')
    plt.text(4.9, 1.9, 'Residual', fontsize=fs2, horizontalalignment='right')





    # plt.show()
    # fig.savefig('fig10.pdf', dpi=300)
    fig.savefig('fig10.png', dpi=300)

    return 0

main()
