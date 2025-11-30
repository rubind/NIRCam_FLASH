from FileRead import readcol, writecol
import matplotlib.pyplot as plt
import numpy as np
import sys


[source_id, x_sw, y_sw, x_lw, y_lw, ra, dec, F150W2_r_small_pix, F150W2_apcorr_mag, F150W2_flux_jy, F150W2_fluxerr_jy, F150W2_mag_ab, F150W2_magerr, F322W2_r_small_pix, F322W2_apcorr_mag, F322W2_flux_jy, F322W2_fluxerr_jy, F322W2_mag_ab, F322W2_magerr, color_f150w2_minus_f322w2] = readcol(sys.argv[1], 'f'*20)

color = F150W2_mag_ab - F322W2_mag_ab

inds = np.where((color < 2)*(F150W2_mag_ab < 27)*(color > -0.75)) #(color < 0.25+(F150W2_mag_ab - 24)*0.25/(28 - 24)))

plt.figure(figsize = (18, 12))
plt.plot(color, F150W2_mag_ab, '.', color = 'b', alpha = 0.1)
plt.plot(color[inds], F150W2_mag_ab[inds], '.', color = 'green')
plt.xlim(-2, 2)
plt.savefig("WD_" + sys.argv[1].split(".")[0] + ".pdf", bbox_inches = 'tight')
plt.close()

[source_id, x_sw, y_sw, x_lw, y_lw, ra, dec, F150W2_r_small_pix, F150W2_apcorr_mag, F150W2_flux_jy, F150W2_fluxerr_jy, F150W2_mag_ab, F150W2_magerr, F322W2_r_small_pix, F322W2_apcorr_mag, F322W2_flux_jy, F322W2_fluxerr_jy, F322W2_mag_ab, F322W2_magerr, color_f150w2_minus_f322w2] = [item[inds] for item in [source_id, x_sw, y_sw, x_lw, y_lw, ra, dec, F150W2_r_small_pix, F150W2_apcorr_mag, F150W2_flux_jy, F150W2_fluxerr_jy, F150W2_mag_ab, F150W2_magerr, F322W2_r_small_pix, F322W2_apcorr_mag, F322W2_flux_jy, F322W2_fluxerr_jy, F322W2_mag_ab, F322W2_magerr, color_f150w2_minus_f322w2]]

writecol("WD_" + sys.argv[1].split(".")[0] + ".txt", [source_id, x_sw, y_sw, x_lw, y_lw, ra, dec, F150W2_r_small_pix, F150W2_apcorr_mag, F150W2_flux_jy, F150W2_fluxerr_jy, F150W2_mag_ab, F150W2_magerr, F322W2_r_small_pix, F322W2_apcorr_mag, F322W2_flux_jy, F322W2_fluxerr_jy, F322W2_mag_ab, F322W2_magerr, color_f150w2_minus_f322w2], headings = ["source_id", "x_sw", "y_sw", "x_lw", "y_lw", "ra", "dec", "F150W2_r_small_pix", "F150W2_apcorr_mag", "F150W2_flux_jy", "F150W2_fluxerr_jy", "F150W2_mag_ab", "F150W2_magerr", "F322W2_r_small_pix", "F322W2_apcorr_mag", "F322W2_flux_jy", "F322W2_fluxerr_jy", "F322W2_mag_ab", "F322W2_magerr", "color_f150w2_minus_f322w2"])
