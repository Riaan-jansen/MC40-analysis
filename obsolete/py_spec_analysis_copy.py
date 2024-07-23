'''Gamma spec analysis tool, .py version.'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
import csv
from spec_utilities import *
# from tkinter import *
# from mpl_interactions import ioff, panhandler, zoom_factory


# importing the data in csv form provided by Cambio
filename = 'Spectra.csv'
file = get_data(filename)
threshold = 15  # threshold difference to reject curve_fit and try again

# print a list of the .spe files to choose from (element is the important bit)
for m in file.spectras:
    for j, char in enumerate(file.spectras):
        print(f'{j} -- {char}', flush=True)

    index = int(input('Input integer of row of chosen file: '))

    specname = file.spectras[index]
    # appends data from the channel columns to form an x-axis
    # appends detector data from relevant column to counts
    for i, line in enumerate(file.fdata):
        if i >= 13:
            file.channels.append(float(line[0]))
            file.counts.append(float(line[index]))

    ltime = float(file.ltimes[index])
    rtime = float(file.rtimes[index])
    lrtime = ltime/rtime

    # basic plot to locate the channel number of the relevant peak
    fig, ax = plotspec(file.channels, file.counts, title=specname)
    plt.show()

    # using scipy curve_fit and spec_utilities file, fits a gaussian
    guess = guess_fit(file.channels, file.counts)
    fit, popt = gauss_fit(file.channels, file.counts, guess)

    # fig, ax = plotspec(file.channels, file.counts, title=specname)
    ax.plot(file.channels, fit, 'r', ls='--', label='Gaussian Curve')
    plt.show()

    # identifying the energy peaks associated with the file isotope
    gammas = g_dict()
    for str in gammas.keys():
        if str in specname:
            peaks = gammas[str]
            element = str
            print(f'List of gamma peaks for {element} isotope:')
    for i, ele in enumerate(peaks):
        print(f'{i} -- {ele[0]}', flush=True)

    # initial calibration
    peak_index = int(input('Select index of peak fitted above:'))
    tuple = peaks[peak_index]
    cal_init = tuple[0]/popt[0]
    energy = [i*cal_init for i in file.channels]
    ppeaks = [k[0] for k in peaks]
    areas = {}
    effs = []

    # setting up plot to add fitted y data for each peak
    fig, axs = plotspec(energy, file.counts, title=specname,
                        xlabel='Energy (keV)')

    for k in peaks:
        # the energy of the peak to calibrate to
        cal_erg = k[0]
        I = k[1]

        guess = guess_fit(energy, file.counts, guess_0=cal_erg)
        fits, popts = gauss_fit(energy, file.counts, guess)

        if np.abs(cal_erg - popts[0]) >= threshold:
            # slicing subroutine
            print(f'{cal_erg} peak is not fitted')
            x_slice, y_slice = slicer(energy, file.counts, guess)
            fit_slice, popts = gauss_fit(x_slice, y_slice, guess)
            axs.plot(x_slice, fit_slice, 'r', ls='--', label=f'{cal_erg} keV peak')
            axs.annotate(f'{cal_erg} keV', (popts[0], popts[1]))
            area_value = simpson(fit_slice, x_slice)
            area_peak = {f'{cal_erg} keV': area_value}
            areas.update(area_peak)
        else:
            # plot individual peaks
            axs.plot(energy, fits, 'r', ls='--', label=f'{cal_erg} keV peak')
            axs.annotate(f'{cal_erg} keV', (popts[0], popts[1]))
            print(f'{cal_erg} keV peak, {popts[0]:.2f} keV : {popts[1]:.2f} counts')
            # numerical integration to calculate area under peak
            area_value = simpson(fits, energy)
            area_peak = {f'{cal_erg} keV': area_value}
            areas.update(area_peak)
        
        eff = efficiency_calc(area_value, rtime, I, element, specname)
        effs.append(eff)
        print(effs)

    efig, eax = plot_eff(effs, ppeaks)

    output = [specname, len(peaks), rtime]

    save = input('Do you want to write to file? y/n: ')
    if save == 'y':
        with open('efficiency.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(output)
    else:
        pass

# displays spectrum and fitted peaks
plt.show()
