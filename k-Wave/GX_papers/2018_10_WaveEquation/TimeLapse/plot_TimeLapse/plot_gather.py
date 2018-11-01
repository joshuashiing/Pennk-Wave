import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def main():

    ############################################################
    # Figure setup

    fx = 4
    fy = 6
    lx = 0.35
    ly = 0.25

    el = (1 - lx * 2) / 2
    et = (1 - ly * 3) / 2
    er = el
    eb = et

    n1 = 450
    n2 = 1201

    x1 = 1635.1
    x2 = 1680.1

    fig = plt.figure(1, figsize=(fx, fy))
    cmap = 'seismic'
    p_lim = 7

    ############################################################
    # Case 1 (time 01)

    data_file = 'case01.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)

    t = data['t_axis']
    d10 = data['d'].transpose()
    d1 = d10[n1 : n2, :]

    ax = fig.add_axes((el, eb + ly * 2, lx, ly))
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Time (ms)')
    ax.imshow(d1, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(x1, x2, t[0, n2-1] * 1e3, t[0, n1] * 1e3), aspect='auto')

    ############################################################
    # Case 2 (time 02)

    data_file = 'case02.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)

    t = data['t_axis']
    d20 = data['d'].transpose()
    d2 = d20[n1: n2, :]

    ax = fig.add_axes((el, eb + ly, lx, ly))
    ax.set_xticks([])
    ax.set_ylabel('Time (ms)')
    ax.imshow(d2, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(x1, x2, t[0, n2 - 1] * 1e3, t[0, n1] * 1e3), aspect='auto')

    ax = fig.add_axes((el + lx, eb + ly, lx, ly))
    ax.set_xticks([])
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.set_ylabel('Time (ms)')
    ax.imshow(d2 - d1, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(x1, x2, t[0, n2 - 1] * 1e3, t[0, n1] * 1e3), aspect='auto')

    ############################################################
    # Case 3 (time 03)

    data_file = 'case03.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)

    t = data['t_axis']
    d30 = data['d'].transpose()
    d3 = d30[n1: n2, :]

    ax = fig.add_axes((el, eb, lx, ly))
    ax.set_xticks([1640, 1660])
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Time (ms)')
    ax.imshow(d3, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(x1, x2, t[0, n2 - 1] * 1e3, t[0, n1] * 1e3), aspect='auto')

    ax = fig.add_axes((el + lx, eb, lx, ly))
    ax.set_xlabel('Distance (m)')
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.set_ylabel('Time (ms)')
    ax.imshow(d3 - d1, cmap=cmap, vmin=-p_lim, vmax=p_lim, extent=(x1, x2, t[0, n2 - 1] * 1e3, t[0, n1] * 1e3), aspect='auto')

    # plt.show()
    fig.savefig('fig_gather.pdf', dpi=300)

    return 0

main()
