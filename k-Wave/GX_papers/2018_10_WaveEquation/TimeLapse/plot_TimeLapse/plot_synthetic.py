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
    lx = 0.75
    ly = 0.8

    el = (1 - lx) *3 / 4
    et = (1 - ly) / 4
    er = 1 - lx - el
    eb = 1 - ly - et

    n1 = 450
    n2 = 1201
    n1 = 550
    n2 = 901

    x1 = 1635.1
    x2 = 1680.1
    x_rec = np.linspace(x1, x2, 151)
    i_rec = [15 * n for n in range(11)]

    # i_rec = np.linspace(0, 150, 11)

    fig = plt.figure(1, figsize=(fx, fy))
    ax = fig.add_axes((el, eb, lx, ly))
    ax.invert_yaxis()
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Distance (m)')

    r = 1

    ############################################################
    # Loading Data

    data_file = 'case01.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    t = data['t_axis'].reshape(-1)
    d1 = data['d']

    data_file = 'case02.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    t = data['t_axis'].reshape(-1)
    d2 = data['d']

    data_file = 'case03.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    t = data['t_axis'].reshape(-1) * 1000
    d3 = data['d']

    data_file = 'case04.mat'
    data_file = os.path.join('..', data_file)
    data = sio.loadmat(data_file)
    t = data['t_axis'].reshape(-1) * 1000
    d4 = data['d']

    x_d = np.expand_dims(x_rec, 1) * np.ones((1, t.shape[0]))
    d1x = -d1 * r + x_d
    d2x = -d2 * r + x_d
    d3x = -d3 * r + x_d
    d4x = -d4 * r + x_d

    for i in i_rec:
        ax.plot(t[n1:n2], d1x[i, n1 : n2], 'k')
        ax.plot(t[n1:n2], d2x[i, n1 : n2], 'b--')
        ax.plot(t[n1:n2], d3x[i, n1 : n2], 'r--')
        # ax.plot(t[n1:n2], d4x[i, n1: n2], 'g--')

    # plt.show()
    fig.savefig('fig_synthetic.pdf', dpi=300)

    return 0

main()
