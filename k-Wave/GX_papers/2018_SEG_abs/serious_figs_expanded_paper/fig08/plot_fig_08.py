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
    t_lim = 0 / 1000
    b_lim = (Nx - 1) * h / 1000

    data = sio.loadmat('BP_cut_model_gas_chimney.mat')
    c = data['velp'] / 1000
    Q = data['Qp']

    ############################################################
    # Plotting data

    cmap1 = 'jet'
    cmap2 = 'viridis'
    # alpha1 = 0.7
    # cmap2 = 'rainbow_r'
    fs = 14

    fig = plt.figure(1, figsize=(7, 6))

    # ax1 = fig.add_axes((.1, .5, .72, .32914))
    ax1 = fig.add_axes((.1, .5, .8, .4))
    ax1.xaxis.set_ticks_position('both')
    ax1.set_xticklabels([])
    ax1.yaxis.set_ticks_position('left')
    plt.ylabel('Distance (km)')
    imgc = ax1.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 1.95, '(a)', fontsize=fs)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad="3%")
    plt.colorbar(imgc, cax=cax, label='Velocity (km/s)')


    ax2 = fig.add_axes((.1, .1, .8, .4))
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('left')
    plt.xlabel('Distance (km)')
    plt.ylabel('Distance (km)')
    imgQ = ax2.imshow(Q, cmap=cmap2, extent=(l_lim, r_lim, b_lim, t_lim))
    plt.text(0.05, 1.95, '(b)', fontsize=fs)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="3%", pad="3%")
    plt.colorbar(imgQ, cax=cax, label=r'$Q$')

    # plt.show()
    fig.savefig('fig08.pdf', dpi=300)

    return 0

main()
