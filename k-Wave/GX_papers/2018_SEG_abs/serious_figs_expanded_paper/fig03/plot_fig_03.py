import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def main():

    ############################################################
    # Loading data

    data1 = sio.loadmat('case01.mat')
    data2 = sio.loadmat('case02.mat')
    data3 = sio.loadmat('case03.mat')
    data4 = sio.loadmat('case04.mat')

    t_axis = data1['t_axis'].reshape(-1)
    dt1 = data1['d2'].reshape(-1)
    da1 = data1['d3'].reshape(-1)
    dt2 = data2['d2'].reshape(-1)
    da2 = data2['d3'].reshape(-1)
    dt3 = data3['d2'].reshape(-1)
    da3 = data3['d3'].reshape(-1)
    dt4 = data4['d2'].reshape(-1)
    da4 = data4['d3'].reshape(-1)

    ##############################################################
    # Plotting

    lw1 = 3
    lw2 = 1
    nf_sd = 30
    color1 = np.array([100., 169., 226.]) / 255.
    color2 = np.array([211., 35., 38.]) / 255.

    plim = 0.055
    fig = plt.figure(1, figsize=(10, 4))
    ax1 = fig.add_axes((.09, .55, .39, .4))
    ax1.xaxis.set_ticks_position('none')
    ax1.set_xticklabels([])
    ax1.set_xlim(0.15, 0.25)
    ax1.set_ylim(-plim, plim)
    plt.ylabel('Pressure')
    plt.plot(t_axis, dt1, color=color1, linewidth=lw1, label='Analytical')
    plt.plot(t_axis, da1, '--', color=color2, linewidth=lw1, label='Numerical')
    plt.text(0.152, -0.050, r'$Q=\infty$ (Acoustic)', fontsize=12)
    ax1.legend(loc='upper left', shadow=True, fontsize=10)


    plim = 0.024
    ax2 = fig.add_axes((.09, .15, .39, .4))
    ax2.set_xlim(0.15, 0.25)
    ax2.set_ylim(-plim, plim)
    plt.xticks(np.linspace(0.15, 0.25, 6))
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure')
    plt.plot(t_axis, dt2, color=color1, linewidth=lw1, label='Analytical')
    plt.plot(t_axis, da2, '--', color=color2, linewidth=lw1, label='Numerical')
    plt.text(0.152, 0.017, r'$Q=100$', fontsize=12)
    ax2.legend(loc='lower left', shadow=True, fontsize=10)

    plim = 0.0055
    ax3 = fig.add_axes((0.52, .55, .39, .4))
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position('right')
    ax3.set_xlim(0.15, 0.25)
    ax3.set_ylim(-plim, plim)
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure')
    plt.plot(t_axis, dt3, color=color1, linewidth=lw1, label='Analytical')
    plt.plot(t_axis, da3, '--', color=color2, linewidth=lw1, label='Numerical')
    plt.text(0.232, -0.005, r'$Q=32$', fontsize=12)
    ax3.legend(loc='upper right', shadow=True, fontsize=10)

    plim = 0.00024
    ax4 = fig.add_axes((.52, .15, .39, .4))
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position('right')
    ax4.set_xlim(0.15, 0.25)
    ax4.set_ylim(-plim, plim)
    plt.xticks(np.linspace(0.15, 0.25, 6))
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure')
    plt.plot(t_axis, dt4, color=color1, linewidth=lw1, label='Analytical')
    plt.plot(t_axis, da4, '--', color=color2, linewidth=lw1, label='Numerical')
    plt.text(0.232, 0.00018, r'$Q=10$', fontsize=12)
    ax4.legend(loc='lower right', shadow=True, fontsize=10)


    # plt.show()
    fig.savefig('fig03.pdf', dpi=300)

    return 0

main()
