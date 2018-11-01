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

    x_src = 0
    y_src = 1658.05
    x0_rec = 30 * np.ones(151)
    y0_rec = np.linspace(1635.1, 1680.1, 151)
    x_rec = 30 * np.ones(11)
    y_rec = np.linspace(1635.1, 1680.1, 11)

    ############################################################
    # Figure setup
    fx = 7
    fy = 7
    lx = 0.25
    et = 0.1

    ly = fx * lx / (Nx - 1) * (Ny - 1) / fy
    eb = 1 - et - ly * 3
    el = (1 - lx * 3) / 2
    er = el

    ec1 = lx * 0.15
    ec2 = ly * 0.33
    lcy = ly * 0.1

    fig = plt.figure(1, figsize=(fx, fy))
    cmap1 = 'jet'
    # cmap2 = 'viridis'
    cmap2 = 'Greens_r'
    cmap3 = 'Reds_r'
    fs = 14

    ms = 12
    mr = 5
    cs = 'k'

    cmin = 2.2
    cmax = 2.8
    Qmin = 0
    Qmax = 50
    dcmin = -400
    dcmax = 0

    ############################################################
    # Case 1 (time 01)

    data_file = 'case01.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c1 = data['vp'] / 1000
    Q1 = data['Q']

    ax = fig.add_axes((el, eb + ly * 2, lx, ly))
    imgQ = ax.imshow(Q1, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xticks([-20, 0, 20])
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_label_position('left')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Distance (m)')

    ax = fig.add_axes((el + lx, eb + ly * 2, lx, ly))
    imgc = ax.imshow(c1, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xticks([-20, 0, 20, 40])
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Distance (m)')

    ############################################################
    # Case 2 (time 02)

    data_file = 'case02.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c2 = data['vp'] / 1000
    Q2 = data['Q']

    ax = fig.add_axes((el, eb + ly, lx, ly))
    imgQ = ax.imshow(Q2, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.set_xticks([])
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_label_position('left')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_ylabel('Distance (m)')

    ax = fig.add_axes((el + lx, eb + ly, lx, ly))
    imgc = ax.imshow(c2, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = fig.add_axes((el + lx * 2, eb + ly, lx, ly))
    imgdc = ax.imshow((c2 - c1) * 1000, cmap=cmap3, vmin=dcmin, vmax=dcmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.set_xticks([])
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_ylabel('Distance (m)')

    ############################################################
    # Case 3 (time 03)

    data_file = 'case03.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c3 = data['vp'] / 1000
    Q3 = data['Q']

    ax = fig.add_axes((el, eb, lx, ly))
    imgQ = ax.imshow(Q3, cmap=cmap2, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_xticks([-20, 0, 20])
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_label_position('left')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Distance (m)')
    # colorbar
    cax = plt.axes((el + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgQ, cax=cax, orientation='horizontal', label='$Q$', ticks=np.linspace(Qmin, Qmax, 3))

    ax = fig.add_axes((el + lx, eb, lx, ly))
    imgc = ax.imshow(c3, cmap=cmap1, vmin=cmin, vmax=cmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_xticks([-20, 0, 20])
    ax.set_yticks([])
    ax.set_xlabel('Distance (m)')
    # colorbar
    cax = plt.axes((el + lx + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgc, cax=cax, orientation='horizontal', label='V (km/s)', ticks=np.linspace(cmin, cmax, 3))

    ax = fig.add_axes((el + lx * 2, eb, lx, ly))
    imgdc = ax.imshow((c3 - c1) * 1000, cmap=cmap3, vmin=dcmin, vmax=dcmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_xticks([-20, 0, 20, 40])
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.set_yticks([1640, 1660, 1680])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Distance (m)')
    # colorbar
    cax = plt.axes((el + lx * 2 + ec1, eb - ec2 - lcy, lx - ec1 * 2, lcy))
    plt.colorbar(imgdc, cax=cax, orientation='horizontal', label='$\Delta$V (m/s)', ticks=np.linspace(dcmin, dcmax, 3))

    # plt.show()
    fig.savefig('fig_models.pdf', dpi=300)

    return 0

main()
