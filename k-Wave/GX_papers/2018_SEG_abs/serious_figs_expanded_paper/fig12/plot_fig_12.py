import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():

    ############################################################
    # Loading data

    data1 = sio.loadmat('case01_k.mat')

    f = data1['f'].reshape(-1)  # Frequency axis
    kt = data1['k'].reshape(-1)  # Theoretical
    ka_111110 = data1['k_c_111110'].reshape(-1)
    ka_111100 = data1['k_c_111100'].reshape(-1)
    ka_111010 =   data1['k_111010'].reshape(-1)
    ka_011100 = data1['k_c_011100'].reshape(-1)
    ka_110100 = data1['k_c_110100'].reshape(-1)


    ############################################################
    # Processing data

    w = 2 * np.pi * f
    cpt = w / np.real(kt) / 1000                # Phase velocity
    alt = -np.imag(kt)                          # attenuation coefficient
    dbt = -20 * np.log10(np.exp(-alt)) * 1000   # alpha in dB/km

    cpa_111110 = w / np.real(ka_111110) / 1000
    ala_111110 = -np.imag(ka_111110)
    dba_111110 = -20 * np.log10(np.exp(-ala_111110)) * 1000

    cpa_111100 = w / np.real(ka_111100) / 1000
    ala_111100 = -np.imag(ka_111100)
    dba_111100 = -20 * np.log10(np.exp(-ala_111100)) * 1000

    cpa_111010 = w / np.real(ka_111010) / 1000
    ala_111010 = -np.imag(ka_111010)
    dba_111010 = -20 * np.log10(np.exp(-ala_111010)) * 1000

    cpa_011100 = w / np.real(ka_011100) / 1000
    ala_011100 = -np.imag(ka_011100)
    dba_011100 = -20 * np.log10(np.exp(-ala_011100)) * 1000

    cpa_110100 = w / np.real(ka_110100) / 1000
    ala_110100 = -np.imag(ka_110100)
    dba_110100 = -20 * np.log10(np.exp(-ala_110100)) * 1000


    ############################################################
    # Plotting data

    lw1 = 2.5
    lw2 = 1
    nf_sd = 30

    color_grid = np.array([220., 220., 220.]) / 255.
    color1 = np.array([100., 169., 226.]) / 255.
    color2 = np.array([211., 35., 38.]) / 255.
    color3 = np.array([221., 198., 186.]) / 255.
    color4 = np.array([207., 158., 118.]) / 255.
    color5 = np.array([146., 109., 90.]) / 255.
    color6 = np.array([84., 62., 54.]) / 255.
    color7 = np.array([112., 112., 112.]) / 255.

    alpha = 0.7



    fig = plt.figure(1, figsize=(10, 4))
    ax1 = fig.add_axes((.08, .4, .4, .55))
    ax1.xaxis.set_ticks_position('none')
    ax1.set_xticklabels([])
    ax1.grid(linestyle='--', linewidth=1, color=color_grid)
    ax1.set_xlim(10, 60)
    plt.plot(f, cpt, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f, cpa_111110, color=color2, linewidth=lw1, label='S111110')
    plt.plot(f, cpa_111100, '--', color=color3, linewidth=lw1, label='S111100')
    plt.plot(f, cpa_111010, '--', color=color4, linewidth=lw1, label='S111010')
    plt.plot(f, cpa_011100, '--', color=color5, linewidth=lw1, label='S011100')
    plt.plot(f, cpa_110100, '--', color=color6, linewidth=lw1, label='S110100')
    plt.ylabel('Velocity (km/s)')
    ax1.legend(loc='lower right', shadow=True, fontsize=8)

    ax2 = fig.add_axes((.08, .15, .4, .25))
    ax2.grid(linestyle='--', linewidth=1, color=color_grid)
    ax2.set_xlim(10, 60)
    f_ds = np.linspace(f[1], f[-1], nf_sd)
    dcpa_111110 = np.interp(f_ds, f, cpa_111110 - cpt)
    dcpa_111100 = np.interp(f_ds, f, cpa_111100 - cpt)
    dcpa_111010 = np.interp(f_ds, f, cpa_111010 - cpt)
    dcpa_011100 = np.interp(f_ds, f, cpa_011100 - cpt)
    dcpa_110100 = np.interp(f_ds, f, cpa_110100 - cpt)
    plt.plot(f_ds, dcpa_111110, 'o--', color=color7, linewidth=lw2, label='S111110', alpha=alpha)
    plt.plot(f_ds, dcpa_111100, '^--', color=color7, linewidth=lw2, label='S111100', alpha=alpha)
    plt.plot(f_ds, dcpa_111010, 's--', color=color7, linewidth=lw2, label='S111010', alpha=alpha)
    plt.plot(f_ds, dcpa_011100, 'D--', color=color7, linewidth=lw2, label='S011100', alpha=alpha)
    plt.plot(f_ds, dcpa_110100, 'X--', color=color7, linewidth=lw2, label='S110100', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Residual (m/s)')
    ax2.legend(ncol=3, loc='lower right', shadow=True, fontsize=8)



    ax3 = fig.add_axes((.52, .4, .4, .55))
    ax3.xaxis.set_ticks_position('none')
    ax3.set_xticklabels([])
    ax3.yaxis.set_ticks_position('right')
    ax3.yaxis.set_label_position('right')
    ax3.grid(linestyle='--', linewidth=1, color=color_grid)
    ax3.set_xlim(10, 60)
    plt.plot(f, dbt, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f, dba_111110, color=color2, linewidth=lw1, label='S111110')
    plt.plot(f, dba_111100, '--', color=color3, linewidth=lw1, label='S111100')
    plt.plot(f, dba_111010, '--', color=color4, linewidth=lw1, label='S111010')
    plt.plot(f, dba_011100, '--', color=color5, linewidth=lw1, label='S011100')
    plt.plot(f, dba_110100, '--', color=color6, linewidth=lw1, label='S110100')
    plt.ylabel(r'$\alpha$ (dB/km)')
    ax3.legend(loc='lower right', shadow=True, fontsize=8)

    tmp_scale = 1
    ax4 = fig.add_axes((.52, .15, .4, .25))
    ax4.yaxis.set_ticks_position('right')
    ax4.yaxis.set_label_position('right')
    ax4.grid(linestyle='--', linewidth=1, color=color_grid)
    ax4.set_xlim(10, 60)
    f_ds = np.linspace(f[1], f[-1], nf_sd)
    ddba_111110 = np.interp(f_ds, f, dba_111110 - dbt) * tmp_scale
    ddba_111100 = np.interp(f_ds, f, dba_111100 - dbt) * tmp_scale
    ddba_111010 = np.interp(f_ds, f, dba_111010 - dbt) * tmp_scale
    ddba_011100 = np.interp(f_ds, f, dba_011100 - dbt) * tmp_scale
    ddba_110100 = np.interp(f_ds, f, dba_110100 - dbt) * tmp_scale
    plt.plot(f_ds, ddba_111110, 'o--', color=color7, linewidth=lw2, label='S111110', alpha=alpha)
    plt.plot(f_ds, ddba_111100, '^--', color=color7, linewidth=lw2, label='S111100', alpha=alpha)
    # plt.plot(f_ds, ddba_111010, 's--', color=color7, linewidth=lw2, label='S111010', alpha=alpha)
    plt.plot(f_ds, ddba_011100, 'D--', color=color7, linewidth=lw2, label='S011100', alpha=alpha)
    plt.plot(f_ds, ddba_110100, 'X--', color=color7, linewidth=lw2, label='S110100', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$\Delta\alpha$ (dB/km)')
    # ax4.legend(ncol=3, loc='lower right', shadow=True, fontsize=8)
    ax4.legend(ncol=2, shadow=True, fontsize=8)


    # plt.show()
    fig.savefig('fig12.pdf', dpi=300)

    return 0

main()
