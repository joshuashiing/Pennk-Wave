import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():

    ############################################################
    # Loading data

    data1a = sio.loadmat('case01_a.mat')
    data1b = sio.loadmat('case01_b.mat')
    data2a = sio.loadmat('case02_a.mat')
    data2b = sio.loadmat('case02_b.mat')
    data3a = sio.loadmat('case03_a.mat')
    data3b = sio.loadmat('case03_b.mat')

    Nx = 200
    Ny = 200
    h = 8.

    l_lim = 0 / 1000
    r_lim = (Ny - 1) * h / 1000
    t_lim = 0 / 1000
    b_lim = (Nx - 1) * h / 1000

    ############################################################
    # Processing data

    p1a = data1a['p']
    p1b = data1b['p']
    p2a = data2a['p']
    p2b = data2b['p']
    p3a = data3a['p']
    p3b = data3b['p']

    ############################################################
    # Plotting data

    ####### Case A

    p_lim = 0.5
    fs = 14
    # cmap = 'coolwarm'
    cmap = 'binary'

    fig = plt.figure(1, figsize=(8, 8))

    ax1 = fig.add_axes((.1, .634, .267, .267))
    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax1.imshow(p1a, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(a)', fontsize=fs)

    ax2 = fig.add_axes((.367, .634, .267, .267))
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_xticks([0.5, 1, 1.5])
    ax2.yaxis.set_ticks_position('none')
    ax2.set_yticks([])
    plt.xlabel('Distance (km)')
    ax2.imshow(p1b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(b)', fontsize=fs)

    ax3 = fig.add_axes((.634, .634, .267, .267))
    ax3.xaxis.set_ticks_position('top')
    ax3.xaxis.set_label_position('top')
    ax3.set_xticks([0.5, 1, 1.5])
    ax3.yaxis.set_ticks_position('right')
    ax3.yaxis.set_label_position('right')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax3.imshow(p1a - p1b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(c)', fontsize=fs)



    ###### Case B

    ax4 = fig.add_axes((.1, .367, .267, .267))
    ax4.xaxis.set_ticks_position('none')
    ax4.set_xticks([])
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')
    plt.ylabel('Distance (km)')
    ax4.imshow(p2a, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(d)', fontsize=fs)

    ax5 = fig.add_axes((.367, .367, .267, .267))
    ax5.xaxis.set_ticks_position('none')
    ax5.set_xticks([])
    ax5.yaxis.set_ticks_position('none')
    ax5.set_yticks([])
    ax5.imshow(p2b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(e)', fontsize=fs)

    ax6 = fig.add_axes((.634, .367, .267, .267))
    ax6.xaxis.set_ticks_position('none')
    ax6.set_xticks([])
    ax6.yaxis.set_ticks_position('right')
    ax6.yaxis.set_label_position('right')
    plt.ylabel('Distance (km)')
    ax6.imshow(p2a - p2b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(f)', fontsize=fs)

    ###### Case C

    ax7 = fig.add_axes((.1, .1, .267, .267))
    ax7.xaxis.set_ticks_position('bottom')
    ax7.xaxis.set_label_position('bottom')
    ax7.yaxis.set_ticks_position('left')
    ax7.yaxis.set_label_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax7.imshow(p3a, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(g)', fontsize=fs)

    ax8 = fig.add_axes((.367, .1, .267, .267))
    ax8.xaxis.set_ticks_position('bottom')
    ax8.xaxis.set_label_position('bottom')
    ax8.set_xticks([0.5, 1, 1.5])
    ax8.yaxis.set_ticks_position('none')
    ax8.set_yticks([])
    plt.xlabel('Distance (km)')
    ax8.imshow(p3b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(h)', fontsize=fs)

    ax9 = fig.add_axes((.634, .1, .267, .267))
    ax9.xaxis.set_ticks_position('bottom')
    ax9.xaxis.set_label_position('bottom')
    ax9.set_xticks([0.5, 1, 1.5])
    ax9.yaxis.set_ticks_position('right')
    ax9.yaxis.set_label_position('right')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax9.imshow(p3a - p3b, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 0.15, '(i)', fontsize=fs)




    # plt.show()
    fig.savefig('fig06.pdf', dpi=300)

    return 0

main()
