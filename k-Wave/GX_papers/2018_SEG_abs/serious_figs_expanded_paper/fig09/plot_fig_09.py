import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():

    ############################################################
    # Loading data

    Nx = 161
    Ny = 398
    h = 12.5
    n_PML = 20


    l_lim = 0 / 1000
    r_lim = (Ny - 1) * h / 1000
    t_lim = 0 / 1000
    b_lim = (Nx - 1) * h / 1000

    data = sio.loadmat('BP_cut_model_gas_chimney.mat')
    data1 = sio.loadmat('case01.mat')
    data2 = sio.loadmat('case02.mat')
    data3 = sio.loadmat('case03.mat')

    c = data['velp'] / 1000
    p1 = data1['p'][n_PML:-n_PML,n_PML:-n_PML]
    p2 = data2['p'][n_PML:-n_PML,n_PML:-n_PML]
    p3 = data3['p'][n_PML:-n_PML,n_PML:-n_PML]

    ############################################################
    # Plotting data

    cmap1 = 'jet'
    cmap2 = 'seismic'
    p_lim = 3
    fs = 14
    alpha=0.75

    fig = plt.figure(1, figsize=(10, 4))
    ax1 = fig.add_axes((.1, .55, .4, .4))
    ax1.xaxis.set_ticks_position('none')
    ax1.set_xticklabels([])
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    plt.ylabel('Distance (km)')
    ax1.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    ax1.imshow(p1, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), alpha=alpha)
    plt.text(0.05, 0.25, '(a)', fontsize=fs)
    plt.text(0.05, 1.95, 'Acoustic', fontsize=fs)

    fig = plt.figure(1, figsize=(10, 4))
    ax2 = fig.add_axes((.5, .55, .4, .4))
    ax2.xaxis.set_ticks_position('none')
    ax2.set_xticklabels([])
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    plt.ylabel('Distance (km)')
    ax2.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    ax2.imshow(p2, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), alpha=alpha)
    plt.text(0.05, 0.25, '(b)', fontsize=fs)
    plt.text(0.05, 1.95, 'Viscoacoustic (New)', fontsize=fs)

    fig = plt.figure(1, figsize=(10, 4))
    ax3 = fig.add_axes((.1, .11, .4, .4))
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax3.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    ax3.imshow(p3, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), alpha=alpha)
    plt.text(0.05, 0.25, '(c)', fontsize=fs)
    plt.text(0.05, 1.95, 'Viscoacoustic (Reference)', fontsize=fs)

    fig = plt.figure(1, figsize=(10, 4))
    ax4 = fig.add_axes((.5, .11, .4, .4))
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('right')
    ax4.yaxis.set_label_position('right')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    ax4.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    ax4.imshow(p2-p3, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), alpha=alpha)
    plt.text(0.05, 0.25, '(d)', fontsize=fs)
    plt.text(0.05, 1.95, 'Residual', fontsize=fs)



    # plt.show()
    fig.savefig('fig09.pdf', dpi=300)


    return 0

main()