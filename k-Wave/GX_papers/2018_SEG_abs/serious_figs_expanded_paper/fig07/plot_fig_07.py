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
    icol = 100
    xaxis = np.arange(Nx) * h / 1000
    xmin = np.min(xaxis)
    xmax = np.max(xaxis)

    ############################################################
    # Processing data

    p1a = data1a['p'][:, icol]
    p1b = data1b['p'][:, icol]
    p2a = data2a['p'][:, icol]
    p2b = data2b['p'][:, icol]
    p3a = data3a['p'][:, icol]
    p3b = data3b['p'][:, icol]

    ############################################################
    # Plotting data

    x_inter = 129.5 * h / 1000
    plim = 1.2
    lw1 = 1.5
    lw2 = 1
    li = 2
    cb = np.array([100., 169., 226.]) / 255.
    ca = np.array([211., 35., 38.]) / 255.
    color_grid = np.array([220., 220., 220.]) / 255.
    cd = '0.3'
    ci = '0.2'
    alpha=0.5
    fs = 14

    ####### Case A

    fig = plt.figure(1, figsize=(8, 5))

    ax1 = fig.add_axes((.15, .65, .8, .25))
    ax1.xaxis.set_ticks_position('both')
    ax1.set_xticklabels([])
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(-plim, plim)
    ax1.grid(linestyle='--', linewidth=1, color=color_grid)
    plt.ylabel('Pressure')
    pltb, = plt.plot(xaxis, p1b, color=cb, linewidth=lw1, label='Reference')
    plta, = plt.plot(xaxis, p1a, '--', color=ca, linewidth=lw1, label='New Scheme')
    pltd, = plt.plot(xaxis, p1a - p1b, ':', color=cd, linewidth=lw2, label='Residual')
    plt.plot([x_inter, x_inter], [-plim, plim], '--', color=ci, linewidth=li, alpha=alpha)
    plt.text(0.03, 0.80, '(a)', fontsize=fs)
    ax1.legend([plta, pltb, pltd], ['New Scheme', 'Reference', 'Residual'],
               ncol=1, loc='lower left', shadow=False, fontsize=7)

    ax2 = fig.add_axes((.15, .40, .8, .25))
    ax2.xaxis.set_ticks_position('both')
    ax2.set_xticklabels([])
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(-plim, plim)
    ax2.grid(linestyle='--', linewidth=1, color=color_grid)
    plt.ylabel('Pressure')
    pltb, = plt.plot(xaxis, p2b, color=cb, linewidth=lw1, label='Reference')
    plta, = plt.plot(xaxis, p2a, '--', color=ca, linewidth=lw1, label='New Scheme')
    pltd, = plt.plot(xaxis, p2a - p2b, ':', color=cd, linewidth=lw2, label='Residual')
    plt.plot([x_inter, x_inter], [-plim, plim], '--', color=ci, linewidth=li, alpha=alpha)
    plt.text(0.03, 0.80, '(b)', fontsize=fs)
    ax2.legend([plta, pltb, pltd], ['New Scheme', 'Reference', 'Residual'],
               ncol=1, loc='lower left', shadow=False, fontsize=7)

    ax3 = fig.add_axes((.15, .15, .8, .25))
    ax3.xaxis.set_ticks_position('both')
    ax3.xaxis.set_label_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(-plim, plim)
    ax3.grid(linestyle='--', linewidth=1, color=color_grid)
    plt.ylabel('Pressure')
    plt.xlabel('Distance (km)')
    pltb, = plt.plot(xaxis, p3b, color=cb, linewidth=lw1, label='Reference')
    plta, = plt.plot(xaxis, p3a, '--', color=ca, linewidth=lw1, label='New Scheme')
    pltd, = plt.plot(xaxis, p3a - p3b, ':', color=cd, linewidth=lw2, label='Residual')
    plt.plot([x_inter, x_inter], [-plim, plim], '--', color=ci, linewidth=li, alpha=alpha)
    plt.text(0.03, 0.80, '(c)', fontsize=fs)
    ax3.legend([plta, pltb, pltd], ['New Scheme', 'Reference', 'Residual'],
               ncol=1, loc='lower left', shadow=False, fontsize=7)


    # plt.show()
    fig.savefig('fig07.pdf', dpi=300)

    return 0

main()
