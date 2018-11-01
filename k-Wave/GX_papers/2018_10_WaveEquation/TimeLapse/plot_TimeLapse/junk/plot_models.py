import numpy as np
import os
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():

    ############################################################
    # Model parameters

    Nx = 467
    Ny = 434
    h = 0.15

    l_lim = -19.95
    r_lim = l_lim + (Ny - 1) * h
    t_lim = 1621.45
    b_lim = t_lim + (Nx - 1) * h

    ############################################################
    # Figure setup
    fx = 8
    fy = 6
    lx = 0.2
    et = 0.05


    ly = fx * lx / (Nx - 1) * (Ny - 1) / fy
    eb = 1 - et - ly * 3
    el = (1 - lx * 4) / 2
    er = el

    ec1 = lx * 0.15
    ec2 = ly * 0.2
    lcy = ly * 0.1


    fig = plt.figure(1, figsize=(fx, fy))
    cmap1 = 'jet'
    cmap2 = 'viridis'
    fs = 14

    cmin = 2.65
    cmax = 2.79
    Qmin = 10
    Qmax = 50
    dcmin = -350 / 1000
    dcmax = 0
    dQmin = -50
    dQmax = -30

    ############################################################
    # Case 1 (time 01)

    data_file = 'case01.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c1 = data['vp'] / 1000
    Q1 = data['Q']

    ax = fig.add_axes((el + lx, eb + ly * 2, lx, ly))
    imgc = ax.imshow(c1, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.set_xticks([])
    ax.set_yticks([])

    ax = fig.add_axes((el + lx * 2, eb + ly * 2, lx, ly))
    imgQ = ax.imshow(Q1, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')

    ############################################################
    # Case 2 (time 02)

    data_file = 'case02.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c2 = data['vp'] / 1000
    Q2 = data['Q']

    ax = fig.add_axes((el + lx, eb + ly, lx, ly))
    imgc = ax.imshow(c2, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')

    ax = fig.add_axes((el + lx * 2, eb + ly, lx, ly))
    imgQ = ax.imshow(Q2, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')

    ax = fig.add_axes((el, eb + ly, lx, ly))
    imgdc = ax.imshow(c2 - c1, cmap=cmap1, vmin=dcmin, vmax=dcmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')

    ax = fig.add_axes((el + lx * 3, eb + ly, lx, ly))
    imgdQ = ax.imshow(Q2 - Q1, cmap=cmap2, vmin=dQmin, vmax=dQmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')

    ############################################################
    # Case 3 (time 03)

    data_file = 'case03.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c3 = data['vp'] / 1000
    Q3 = data['Q']

    ax = fig.add_axes((el + lx, eb, lx, ly))
    imgc = ax.imshow(c3, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    cax = plt.axes((el + lx + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgc, cax=cax, orientation='horizontal', label='V (km/s)', ticks=np.linspace(cmin, cmax, 3))


    ax = fig.add_axes((el + lx * 2, eb, lx, ly))
    imgQ = ax.imshow(Q3, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    cax = plt.axes((el + lx * 2 + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgQ, cax=cax, orientation='horizontal', label='$Q$', ticks=np.linspace(Qmin, Qmax, 3))

    ax = fig.add_axes((el, eb, lx, ly))
    imgdc = ax.imshow((c3 - c1) * 1000, cmap=cmap1, vmin=dcmin*1000, vmax=dcmax*1000, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    cax = plt.axes((el + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgdc, cax=cax, orientation='horizontal', label='$\Delta$V (m/s)', ticks=np.linspace(dcmin * 1000, dcmax * 1000, 3))

    ax = fig.add_axes((el + lx * 3, eb, lx, ly))
    imgdQ = ax.imshow(Q3 - Q1, cmap=cmap2, vmin=dQmin, vmax=dQmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    cax = plt.axes((el + lx * 3 + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgdQ, cax=cax, orientation='horizontal', label='$\Delta$$Q$', ticks=np.linspace(dQmin, dQmax, 3))

    plt.show()


    # # ax1 = fig.add_axes((.1, .5, .72, .32914))
    # ax1 = fig.add_axes((.1, .5, .8, .4))
    # ax1.xaxis.set_ticks_position('both')
    # ax1.set_xticklabels([])
    # ax1.yaxis.set_ticks_position('left')
    # plt.ylabel('Distance (km)')
    # imgc = ax1.imshow(c, cmap=cmap1, extent=(l_lim, r_lim, b_lim, t_lim))
    # plt.text(0.05, 1.95, 'Velocity Model', fontsize=fs)
    # divider = make_axes_locatable(ax1)
    # cax = divider.append_axes("right", size="3%", pad="3%")
    # plt.colorbar(imgc, cax=cax, label='Velocity (km/s)')
    #
    #
    # ax2 = fig.add_axes((.1, .1, .8, .4))
    # ax2.xaxis.set_ticks_position('both')
    # ax2.yaxis.set_ticks_position('left')
    # plt.xlabel('Distance (km)')
    # plt.ylabel('Distance (km)')
    # imgQ = ax2.imshow(Q, cmap=cmap2, extent=(l_lim, r_lim, b_lim, t_lim))
    # plt.text(0.05, 1.95, '$Q$ Model', fontsize=fs)
    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes("right", size="3%", pad="3%")
    # plt.colorbar(imgQ, cax=cax, label=r'$Q$')
    #
    # # plt.show()
    # fig.savefig('fig08_poster.pdf', dpi=300)

    return 0

main()
