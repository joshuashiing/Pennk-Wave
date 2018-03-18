import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def main():
    # Loading data
    data1 = sio.loadmat('case01_k.mat')
    data2 = sio.loadmat('case02_k.mat')
    data3 = sio.loadmat('case03_k.mat')
    data4 = sio.loadmat('case04_k.mat')
    data5 = sio.loadmat('case05_k.mat')
    data6 = sio.loadmat('case06_k.mat')

    # Case A (10 - 60 Hz)
    f_a = data1['f'].reshape(-1)                # Frequency axis
    kt_a_1 = data1['k'].reshape(-1)             # Theoretical
    ka_a_1 = data1['k_c_111110'].reshape(-1)    # Approximate
    kt_a_2 = data2['k'].reshape(-1)
    ka_a_2 = data2['k_c_111110'].reshape(-1)
    kt_a_3 = data3['k'].reshape(-1)
    ka_a_3 = data3['k_c_111110'].reshape(-1)

    # Case B (10 - 200 Hz)
    f_b = data4['f'].reshape(-1)
    kt_b_1 = data4['k'].reshape(-1)
    ka_b_1 = data4['k_c_111110'].reshape(-1)
    kt_b_2 = data5['k'].reshape(-1)
    ka_b_2 = data5['k_c_111110'].reshape(-1)
    kt_b_3 = data6['k'].reshape(-1)
    ka_b_3 = data6['k_c_111110'].reshape(-1)

    # Process Loaded data
    w_a = 2 * np.pi * f_a
    cpt_a_1 = w_a / np.real(kt_a_1)                     # Phase velocity
    alt_a_1 = -np.imag(kt_a_1)                          # attenuation coefficient
    dbt_a_1 = -20 * np.log10(np.exp(-alt_a_1)) * 1000   # alpha in dB/km
    cpt_a_2 = w_a / np.real(kt_a_2)
    alt_a_2 = -np.imag(kt_a_2)
    dbt_a_2 = -20 * np.log10(np.exp(-alt_a_2)) * 1000
    cpt_a_3 = w_a / np.real(kt_a_3)
    alt_a_3 = -np.imag(kt_a_3)
    dbt_a_3 = -20 * np.log10(np.exp(-alt_a_3)) * 1000

    cpa_a_1 = w_a / np.real(ka_a_1)
    ala_a_1 = -np.imag(ka_a_1)
    dba_a_1 = -20 * np.log10(np.exp(-ala_a_1)) * 1000
    cpa_a_2 = w_a / np.real(ka_a_2)
    ala_a_2 = -np.imag(ka_a_2)
    dba_a_2 = -20 * np.log10(np.exp(-ala_a_2)) * 1000
    cpa_a_3 = w_a / np.real(ka_a_3)
    ala_a_3 = -np.imag(ka_a_3)
    dba_a_3 = -20 * np.log10(np.exp(-ala_a_3)) * 1000

    w_b = 2 * np.pi * f_b
    cpt_b_1 = w_b / np.real(kt_b_1)  # Phase velocity
    alt_b_1 = -np.imag(kt_b_1)  # attenuation coefficient
    dbt_b_1 = -20 * np.log10(np.exp(-alt_b_1)) * 1000
    cpt_b_2 = w_b / np.real(kt_b_2)
    alt_b_2 = -np.imag(kt_b_2)
    dbt_b_2 = -20 * np.log10(np.exp(-alt_b_2)) * 1000
    cpt_b_3 = w_b / np.real(kt_b_3)
    alt_b_3 = -np.imag(kt_b_3)
    dbt_b_3 = -20 * np.log10(np.exp(-alt_b_3)) * 1000

    cpa_b_1 = w_b / np.real(ka_b_1)
    ala_b_1 = -np.imag(ka_b_1)
    dba_b_1 = -20 * np.log10(np.exp(-ala_b_1)) * 1000
    cpa_b_2 = w_b / np.real(ka_b_2)
    ala_b_2 = -np.imag(ka_b_2)
    dba_b_2 = -20 * np.log10(np.exp(-ala_b_2)) * 1000
    cpa_b_3 = w_b / np.real(ka_b_3)
    ala_b_3 = -np.imag(ka_b_3)
    dba_b_3 = -20 * np.log10(np.exp(-ala_b_3)) * 1000

    ##################################################################################
    # Plot section

    lw1 = 3
    lw2 = 1
    nf_sd = 30

    color1 = np.array([100., 169., 226.]) / 255.
    color2 = np.array([211., 35., 38.]) / 255.
    color3 = np.array([112., 112., 112.]) / 255.
    color_grid = np.array([220., 220., 220.]) / 255.
    alpha = 0.7


    ########### Case A | Phase velocity

    fig = plt.figure(1, figsize=(10, 4))
    ax1 = fig.add_axes((.08, .4, .4, .55))
    ax1.xaxis.set_ticks_position('none')
    ax1.set_xticklabels([])
    ax1.grid(linestyle='--', linewidth=1, color=color_grid)
    ax1.set_xlim(10, 60)
    plt.plot(f_a, cpt_a_1 / 1000, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f_a, cpa_a_1 / 1000, '--', color=color2, linewidth=lw1, label='Approximation')
    plt.plot(f_a, cpt_a_2 / 1000, color=color1, linewidth=lw1)
    plt.plot(f_a, cpa_a_2 / 1000, '--', color=color2, linewidth=lw1)
    plt.plot(f_a, cpt_a_3 / 1000, color=color1, linewidth=lw1)
    plt.plot(f_a, cpa_a_3 / 1000, '--', color=color2, linewidth=lw1)
    plt.ylabel('Velocity (km/s)')
    legend = ax1.legend(loc='lower right', shadow=False)
    plt.text(16, 2.02, r'$Q=10$', fontsize=11)
    plt.text(14, 2.103, r'$Q=32$', fontsize=11)
    plt.text(10.1, 2.133, r'$Q=100$', fontsize=11)


    ax2 = fig.add_axes((.08, .15, .4, .25))
    ax2.grid(linestyle='--', linewidth=1, color=color_grid)
    ax2.set_xlim(10, 60)
    f_a_ds = np.linspace(f_a[1], f_a[-1], nf_sd)
    dcp_a_1 = np.interp(f_a_ds, f_a, cpa_a_1 - cpt_a_1)
    dcp_a_2 = np.interp(f_a_ds, f_a, cpa_a_2 - cpt_a_2)
    dcp_a_3 = np.interp(f_a_ds, f_a, cpa_a_3 - cpt_a_3)
    plt.plot(f_a_ds, dcp_a_1, 'o--', color=color3, linewidth=lw2, label=r'$Q=10$', alpha=alpha)
    plt.plot(f_a_ds, dcp_a_2, '^--', color=color3, linewidth=lw2, label=r'$Q=32$', alpha=alpha)
    plt.plot(f_a_ds, dcp_a_3, 's--', color=color3, linewidth=lw2, label=r'$Q=100$', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Residual (m/s)')
    legend = ax2.legend(ncol=3, loc='lower right', shadow=False, fontsize=8)



    ########### Case A | Attenuation coefficient
    ax3 = fig.add_axes((.52, .4, .4, .55))
    ax3.grid(linestyle='--', linewidth=1, color=color_grid)
    ax3.xaxis.set_ticks_position('none')
    ax3.yaxis.tick_right()
    ax3.set_xticklabels([])
    ax3.set_xlim(10, 60)
    ax3.yaxis.set_label_position('right')
    plt.plot(f_a, dbt_a_1, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f_a, dba_a_1, '--', color=color2, linewidth=lw1, label='Approximation')
    plt.plot(f_a, dbt_a_2, color=color1, linewidth=lw1)
    plt.plot(f_a, dba_a_2, '--', color=color2, linewidth=lw1)
    plt.plot(f_a, dbt_a_3, color=color1, linewidth=lw1)
    plt.plot(f_a, dba_a_3, '--', color=color2, linewidth=lw1)
    plt.ylabel(r'$\alpha$ (dB/km)')
    legend = ax3.legend(loc='upper left', shadow=False)
    plt.text(48, 55, r'$Q=10$', fontsize=11)
    plt.text(45, 23, r'$Q=32$', fontsize=11)
    plt.text(50, 0, r'$Q=100$', fontsize=11)


    tmp_scale = 1e5
    ax4 = fig.add_axes((.52, .15, .4, .25))
    ax4.grid(linestyle='--', linewidth=1, color=color_grid)
    ax4.yaxis.tick_right()
    ax4.set_xlim(10, 60)
    ax4.set_ylim(-0.1, 1.7)
    ax4.yaxis.set_label_position('right')
    f_a_ds = np.linspace(f_a[1], f_a[-1], nf_sd)
    dal_a_1 = np.interp(f_a_ds, f_a, ala_a_1 - alt_a_1) * tmp_scale
    dal_a_2 = np.interp(f_a_ds, f_a, ala_a_2 - alt_a_2) * tmp_scale
    dal_a_3 = np.interp(f_a_ds, f_a, ala_a_3 - alt_a_3) * tmp_scale
    plt.plot(f_a_ds, dal_a_1, 'o--', color=color3, linewidth=lw2, label=r'$Q=10$', alpha=alpha)
    plt.plot(f_a_ds, dal_a_2, '^--', color=color3, linewidth=lw2, label=r'$Q=32$', alpha=alpha)
    plt.plot(f_a_ds, dal_a_3, 's--', color=color3, linewidth=lw2, label=r'$Q=100$', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$\Delta\alpha$ ($10^{-5}$ dB/km)')
    legend = ax4.legend(ncol=3, loc='upper right', shadow=False, fontsize=8)

    fig.savefig('fig02a.pdf', dpi=300)




    ########### Case B | Phase velocity

    fig = plt.figure(2, figsize=(10, 4))
    ax1 = fig.add_axes((.08, .4, .4, .55))
    ax1.xaxis.set_ticks_position('none')
    ax1.set_xticklabels([])
    ax1.grid(linestyle='--', linewidth=1, color=color_grid)
    ax1.set_xlim(10, 200)
    plt.plot(f_b, cpt_b_1 / 1000, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f_b, cpa_b_1 / 1000, '--', color=color2, linewidth=lw1, label='Approximation')
    plt.plot(f_b, cpt_b_2 / 1000, color=color1, linewidth=lw1)
    plt.plot(f_b, cpa_b_2 / 1000, '--', color=color2, linewidth=lw1)
    plt.plot(f_b, cpt_b_3 / 1000, color=color1, linewidth=lw1)
    plt.plot(f_b, cpa_b_3 / 1000, '--', color=color2, linewidth=lw1)
    plt.ylabel('Velocity (km/s)')
    legend = ax1.legend(loc='lower right', shadow=False)
    plt.text(30, 2.02, r'$Q=10$', fontsize=11)
    plt.text(14, 2.10, r'$Q=32$', fontsize=11)
    plt.text(10.2, 2.17, r'$Q=100$', fontsize=11)


    ax2 = fig.add_axes((.08, .15, .4, .25))
    ax2.grid(linestyle='--', linewidth=1, color=color_grid)
    ax2.set_xlim(10, 200)
    f_b_ds = np.linspace(f_b[1], f_b[-1], nf_sd)
    dcp_b_1 = np.interp(f_b_ds, f_b, cpa_b_1 - cpt_b_1)
    dcp_b_2 = np.interp(f_b_ds, f_b, cpa_b_2 - cpt_b_2)
    dcp_b_3 = np.interp(f_b_ds, f_b, cpa_b_3 - cpt_b_3)
    plt.plot(f_b_ds, dcp_b_1, 'o--', color=color3, linewidth=lw2, label=r'$Q=10$', alpha=alpha)
    plt.plot(f_b_ds, dcp_b_2, '^--', color=color3, linewidth=lw2, label=r'$Q=32$', alpha=alpha)
    plt.plot(f_b_ds, dcp_b_3, 's--', color=color3, linewidth=lw2, label=r'$Q=100$', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Residual (m/s)')
    legend = ax2.legend(ncol=3, loc='lower right', shadow=False, fontsize=8)




    ########### Case B | Attenuation coefficient
    ax3 = fig.add_axes((.52, .4, .4, .55))
    ax3.grid(linestyle='--', linewidth=1, color=color_grid)
    ax3.xaxis.set_ticks_position('none')
    ax3.yaxis.tick_right()
    ax3.set_xticklabels([])
    ax3.set_xlim(10, 200)
    ax3.yaxis.set_label_position('right')
    plt.plot(f_b, dbt_b_1, color=color1, linewidth=lw1, label='Theoretical')
    plt.plot(f_b, dba_b_1, '--',color=color2, linewidth=lw1, label='Approximation')
    plt.plot(f_b, dbt_b_2, color=color1, linewidth=lw1)
    plt.plot(f_b, dba_b_2, '--', color=color2, linewidth=lw1)
    plt.plot(f_b, dbt_b_3, color=color1, linewidth=lw1)
    plt.plot(f_b, dba_b_3, '--', color=color2, linewidth=lw1)
    plt.ylabel(r'$\alpha$ (dB/km)')
    legend = ax3.legend(loc='upper left', shadow=False)
    plt.text(160, 180, r'$Q=10$', fontsize=11)
    plt.text(145, 70, r'$Q=32$', fontsize=11)
    plt.text(155, 0, r'$Q=100$', fontsize=11)


    tmp_scale = 1e5
    ax4 = fig.add_axes((.52, .15, .4, .25))
    ax4.grid(linestyle='--', linewidth=1, color=color_grid)
    ax4.yaxis.tick_right()
    ax4.set_xlim(10, 200)
    ax4.yaxis.set_label_position('right')
    f_b_ds = np.linspace(f_b[1], f_b[-1], nf_sd)
    dal_b_1 = np.interp(f_b_ds, f_b, ala_b_1 - alt_b_1) * tmp_scale
    dal_b_2 = np.interp(f_b_ds, f_b, ala_b_2 - alt_b_2) * tmp_scale
    dal_b_3 = np.interp(f_b_ds, f_b, ala_b_3 - alt_b_3) * tmp_scale
    plt.plot(f_b_ds, dal_b_1, 'o--', color=color3, linewidth=lw2, label=r'$Q=10$', alpha=alpha)
    plt.plot(f_b_ds, dal_b_2, '^--', color=color3, linewidth=lw2, label=r'$Q=32$', alpha=alpha)
    plt.plot(f_b_ds, dal_b_3, 's--', color=color3, linewidth=lw2, label=r'$Q=100$', alpha=alpha)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$\Delta\alpha$ ($10^{-5}$ dB/km)')
    legend = ax4.legend(ncol=3, loc='upper right', shadow=False, fontsize=8)

    fig.savefig('fig02b.pdf', dpi=300)
    # plt.show()

    return 0

main()
