import numpy as np
import os
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.colors import LinearSegmentedColormap


def main():

    ############################################################
    # Model parameters

    Nx = 467
    Ny = 434
    h = 0.15
    n_PML = 20

    l_lim = -19.95
    r_lim = l_lim + (Ny - 1) * h
    t_lim = 1621.45
    b_lim = t_lim + (Nx - 1) * h

    x_src = 0
    # y_src = 1658.05
    y_src = 1650.10
    x_rec = 30 * np.ones(11)
    y_rec = np.linspace(1635.1, 1680.1, 11)

    ############################################################
    # Figure setup

    fx = 3
    fy = 6
    lx = 0.6
    et = 0.05

    ly = fx * lx / (Nx - 1) * (Ny - 1) / fy
    eb = 1 - et - ly * 3
    el = (1 - lx) * 3/4
    er = el


    Qmin = 0
    Qmax = 50

    p_lim = 2.0
    alpha1 = 0.3
    alpha2 = 0.4
    fs = 14

    ms = 12
    mr = 5
    cs = 'k'

    cmap1 = 'Greens_r'

    rr1 = 4 / 256
    gg1 = 99 / 256
    bb1 = 128 / 256
    rr2 = 255 / 256
    gg2 = 176 / 256
    bb2 = 59 / 256
    cdict_p = {'red': ((0.0, rr1, rr1),
                       (0.5, 1, 1),
                       (1.0, rr2, rr2)),
               'green': ((0.0, gg1, gg1),
                         (0.5, 1, 1),
                         (1.0, gg2, gg2)),
               'blue': ((0.0, bb1, bb1),
                        (0.5, 1, 1),
                        (1.0, bb2, bb2))}
    cmap2 = LinearSegmentedColormap('cmap_p', cdict_p, 256)

    fig = plt.figure(1, figsize=(fx, fy))

    ############################################################
    # Case 1 (time 01)

    data_file = 'case01.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c1 = data['vp'] / 1000
    Q1 = data['Q']
    p1 = data['p'][n_PML:-n_PML, n_PML:-n_PML]

    ax = fig.add_axes((el, eb + ly * 2, lx, ly))
    ax.set_xticks([])
    ax.set_ylabel('Distance (m)')
    ax.imshow(p1, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms, alpha=alpha2)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr, alpha=alpha2)
    ax.imshow(Q1, cmap=cmap1, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto', alpha=alpha1)

    ############################################################
    # Case 2 (time 02)

    data_file = 'case02.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c2 = data['vp'] / 1000
    Q2 = data['Q']
    p2 = data['p'][n_PML:-n_PML, n_PML:-n_PML]

    ax = fig.add_axes((el, eb + ly, lx, ly))
    ax.set_xticks([])
    ax.set_ylabel('Distance (m)')
    ax.imshow(p2, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms, alpha=alpha2)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr, alpha=alpha2)
    ax.imshow(Q2, cmap=cmap1, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto', alpha=alpha1)

    ############################################################
    # Case 3 (time 03)

    data_file = 'case03.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    c3 = data['vp'] / 1000
    Q3 = data['Q']
    p3 = data['p'][n_PML:-n_PML, n_PML:-n_PML]

    ax = fig.add_axes((el, eb, lx, ly))
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Distance (m)')
    ax.imshow(p3, cmap=cmap2, vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto')
    plt.plot(x_src, y_src, '*', color=cs, markersize=ms, alpha=alpha2)
    plt.plot(x_rec, y_rec, 'v', color=cs, markersize=mr, alpha=alpha2)
    ax.imshow(Q3, cmap=cmap1, vmin=Qmin, vmax=Qmax, extent=(l_lim, r_lim, b_lim, t_lim), aspect='auto', alpha=alpha1)

    # plt.show()
    fig.savefig('fig_snapshot.pdf', dpi=300)

    return 0

main()
