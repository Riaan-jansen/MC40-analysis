'''efficiency curve for detector'''
import spec_utilities as su
from spec_utilities import dicts
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.optimize import curve_fit
import argparse
from matplotlib.widgets import Slider, Button
import numpy as np


def peak_lookup(filename):
    '''given spectra name, looks up g_dict and returns peaks and element'''
    gammas = dicts.g_dict()
    for str in gammas.keys():
        if str in filename:
            peaks = gammas[str]
            element = str
    return peaks, element


def write2file(ifile, ofile='efficiency.csv'):
    x, y = main(ifile)
    if len(x) != len(y):
        raise ValueError('x and y not the same length for file writing')
    with open(ofile, 'w', newline='') as file:
        file.write('Energy (keV), Efficiency (abs)\n')
        i = 0
        while i < len(x):
            file.write(f"{x[i]:.2f},{y[i]:.3f},{0.01/(i+1):.3f}\n")
            i = i + 1


def calibration(filename, file):
    '''calibrates the channels to energy, using short detector Eu source
    with exception for spectra binned in very different channels.
    in: spectra name, file class. out: energy: list'''
    cals = dicts.cal_dict()
    for str in cals.keys():
        if str in filename:
            channel_guess = cals[str]
            channels, counts = su.load_data(filename, file)
            break
        elif str not in filename:
            calfile = su.get_data('Spectra.csv')
            channels, counts = su.load_data('82mm_New_152Eu.Spe', calfile)
            filename = '82mm_New_152Eu.Spe'
            channel_guess = (4090, 7)

    # identifying the energy peaks associated with the file isotope
    peaks, ele = peak_lookup(filename)

    # scipy curve_fit and spec_utilities fits a gaussian
    guess = su.guess_fit(channels, counts, guess_0=channel_guess[0])
    fit, popt = su.gauss_fit(channels, counts, guess)

    # calibration using known energy from gammas dict
    peak_index = channel_guess[1]
    tuple = peaks[peak_index]
    cal_init = tuple[0]/popt[0]
    energy = [i*cal_init for i in channels]
    return energy


def main(input, show_plots=False, threshold=10):
    '''executes data reading, calibrates and plots efficiency curve'''
    # loading the data
    file = su.get_data(input)
    forbidden = ['82mm_Al.Spe', 'Background_3-6Nov.Spe', '23cm_152Eu.Spe']
    tuples_cm = []
    tuples_mm = []

    for j, specname in enumerate(file.spectras[1:]):
        if specname in forbidden:
            pass
        else:
            channel, count = su.load_data(specname, file)
            energy = calibration(specname, file)
            peaks, element = peak_lookup(specname)
            if show_plots is True:
                # setting up plot to add fitted y data for each peak
                fig, axs = su.plotspec(energy, count, title=specname,
                                       xlabel='Energy (keV)')
            for k in peaks:
                p_erg = k[0]
                intensity = k[1]

                guess = su.guess_fit(energy, count, guess_0=p_erg)
                fit, popt = su.gauss_fit(energy, count, guess)

                if abs(p_erg - popt[0]) >= threshold:
                    # slicing routine in event peaks too close
                    x_slice, y_slice = su.slicer(energy, count, guess)
                    fit_slice, popts = su.gauss_fit(x_slice, y_slice, guess)
                    area = simpson(fit_slice, x_slice)
                    if show_plots is True:
                        print(f'{p_erg} peak not fitted. Trying refit.')
                        axs.plot(x_slice, fit_slice, 'r', ls='--')
                        axs.annotate(f'{p_erg} keV', (popts[0], popts[1]))

                else:
                    # numerical integration to calculate area under peak
                    area = simpson(fit, energy)

                    if show_plots is True:
                        print(f'{p_erg} peak, {popt[0]:.2f} keV fit')
                        # plot individual peaks
                        axs.plot(energy, fit, 'r', ls='--')
                        axs.annotate(f'{p_erg} keV', (popt[0], popt[1]))

                params = area, file.rtimes[j], intensity, element, specname
                eff = su.efficiency_calc(*params)

                units = dicts.unit_dict()
                for ele in units:
                    if ele in specname:
                        unit_factor = units[ele]
                        if unit_factor == 1:
                            tuples_cm.append((p_erg, eff))
                        else:
                            tuples_mm.append((p_erg, eff))

        if show_plots is True:
            plt.show()

    energy, efficiency = su.unpacker(tuples_cm)
    # su.eff_curve(energy, efficiency)
    # energy_s, efficiency_s = su.unpacker(tuples_mm)
    # su.eff_curve(energy_s, efficiency_s, d='short')

    return energy, efficiency


def fitting(function, guess=None):
    x, y = main('Spectra.csv')
    x = np.asarray(x)
    opt, cov = curve_fit(function, x, y, p0=guess,
                         sigma=0.005*np.ones_like(y))
    print('popts:', opt)
    fit = function(x, *opt)
    plt.scatter(x, y)
    plt.plot(x, fit, color='red')
    plt.show()
    return opt


def sliders():
    def update(val):
        a = amp.val
        b = lamb.val
        c = centr.val
        frame.set_ydata(su.functions.expo(x, a, b, c))

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.35)

    x, y = main('Spectra.csv')
    plt.plot(x, y)
    plt.show()
    x = np.asarray(x)

    guess_e = [0.0009, -0.0041, 200]
    popt = fitting(su.functions.expo, guess_e)
    mfit = su.functions.expo(x, *popt)
    frame, = plt.plot(x, y)

    axamp = plt.axes([0.25, 0.15, 0.65, 0.03])
    axB = plt.axes([0.25, 0.1, 0.65, 0.03])
    axC = plt.axes([0.25, 0.05, 0.65, 0.03])

    amp = Slider(axamp, 'A', 0.00001, 1.0, valinit=0.0009, valstep=0.0001)
    lamb = Slider(axB, 'B', -0.1, -0.001, valinit=-0.01, valstep=0.0001)
    centr = Slider(axC, 'C', -500.0, 500.0, valinit=5.0, valstep=10.0)

    amp.on_changed(update)
    lamb.on_changed(update)
    centr.on_changed(update)
    plt.show()


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="plots efficiency curve")
#     parser.add_argument("input", help="path to Cambio input")
#     parser.add_argument("show_plots", help="show spectra plots")
#     parser.add_argument("write2file", help="boolean, write to file y/n")
#     parser.add_argument("output", help="""to write to file,
#                         (provided input) output filepath""")
#     args = parser.parse_args()

#     main(args.input, args.show_plots)
#     if args.write2file is True:
#         write2file(args.input, args.output)
