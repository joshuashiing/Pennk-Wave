import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():
    ############################################################
    # Loading data

    data1 = sio.loadmat('testS_case01.mat')
    data2 = sio.loadmat('testS_case02.mat')
    data3 = sio.loadmat('testS_case03.mat')
    data4 = sio.loadmat('testS_case04.mat')

    h = 1.

    ############################################################
    # Processing data

    p1 = data1['p']
    p_FQ = np.zeros(p1.shape)
    [Nx, Ny] = p1.shape
    NQx = int(np.floor(Nx / 2))
    NQy = int(np.floor(Ny / 2))
    p_FQ[: NQx, : NQy] = p1[: NQx, : NQy]

    p2 = data2['p']
    p_FQ[: NQx, NQy:] = p2[: NQx, NQy:]

    p3 = data3['p']
    p_FQ[NQx:, :NQy] = p3[NQx:, :NQy]

    p4 = data4['p']
    p_FQ[NQx:, NQy:] = p4[NQx:, NQy:]

    ############################################################
    # Plotting data

    lw = 2

    fig = plt.figure(1, figsize=(6, 6))
    ax = fig.add_axes((.1, .1, .8, .8))

    p_lim = np.max(np.abs(p_FQ))
    l_lim = 0
    r_lim = h * (Ny - 1)
    t_lim = 0
    b_lim = h * (Nx - 1)

    x_mid = h * (NQx - 1)
    y_mid = h * (NQy - 1)

    # ax.imshow(p_FQ, cmap='coolwarm', vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))
    ax.imshow(p_FQ, cmap='binary', vmin=-p_lim, vmax=p_lim, extent=(l_lim, r_lim, b_lim, t_lim))

    plt.plot([l_lim, r_lim], [x_mid, x_mid], '--', linewidth=lw, color='0.2')
    plt.plot([y_mid, y_mid], [t_lim, b_lim], '--', linewidth=lw, color='0.2')

    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')
    plt.text(10, 30, r'(a) Acoustic', fontsize=12)
    plt.text(330, 30, r'(b) Loss-dominated', fontsize=12)
    plt.text(10, 500, r'(c) Dispersion-dominated', fontsize=12)
    plt.text(350, 500, r'(d) Viscoacoustic', fontsize=12)

    # plt.show()
    # fig.savefig('fig05.pdf', dpi=300)
    fig.savefig('fig05.png', dpi=300)

    return 0


main()
