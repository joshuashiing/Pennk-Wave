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
    irec = 280

    l_lim = 0 / 1000
    r_lim = (Ny - 1) * h / 1000
    t_lim = 0.
    b_lim = 2.

    data1 = sio.loadmat('case01.mat')
    data2 = sio.loadmat('case02.mat')
    data3 = sio.loadmat('case03.mat')

    t_axis = data1['t_axis'].reshape(-1)
    d1 = data1['d'][irec, :]
    d2 = data2['d'][irec, :]
    d3 = data3['d'][irec, :]

    ############################################################
    # Plotting data

    cmap = 'seismic'
    p_lim1 = 1.8
    p_lim2 = 0.32
    fs1 = 14
    fs2 = 10
    x1 = 0.85
    x2 = 2.0
    lw1 = 2
    lw2 = 1

    color1 = np.array([100., 169., 226.]) / 255.
    color2 = np.array([211., 35., 38.]) / 255.
    color3 = np.array([112., 112., 112.]) / 255.
    color_grid = np.array([220., 220., 220.]) / 255.
    alpha = 0.7

    fig = plt.figure(1, figsize=(7, 6))
    ax1 = fig.add_axes((.1, .53, .85, .4))
    ax1.spines['right'].set_color(None)
    ax1.spines['top'].set_color(None)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xticklabels([])
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    ax1.set_xlim(x1, x2)
    ax1.set_ylim(-p_lim1, p_lim1)
    plt.ylabel('Pressure')

    plt.plot(t_axis, d1, color=color1, linewidth=lw1, label='Acoustic')
    plt.plot(t_axis, d2, '--', color=color2, linewidth=lw1, label='Viscoacoustic (New)')
    plt.text(.86, 1.4, '(a)', fontsize=fs1)
    ax1.legend(loc='lower right', shadow=False, fontsize=fs2)



    ax2 = fig.add_axes((.1, .1, .85, .4))
    ax2.spines['right'].set_color(None)
    ax2.spines['top'].set_color(None)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.set_xlim(x1, x2)
    ax2.set_ylim(-p_lim2, p_lim2)
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure')

    plt.plot(t_axis, d3, color=color1,  linewidth=lw1, label='Viscoacoustic (Reference)')
    plt.plot(t_axis, d2, '--', color=color2, linewidth=lw1, label='Viscoacoustic (New)')
    plt.plot(t_axis, d2 - d3, ':', color=color3, linewidth=lw2, label='Residual')
    plt.text(.86, 0.26, '(b)', fontsize=fs1)
    ax2.legend(loc='lower right', shadow=False, fontsize=fs2)



    # plt.show()
    fig.savefig('fig11.pdf', dpi=300)

    return 0

main()
