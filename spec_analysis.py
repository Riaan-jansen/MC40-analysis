'''Gamma spec analysis tool, .py version.'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks, find_peaks_cwt
import spec_utilities as su
from spec_utilities import dicts
from matplotlib.widgets import Slider


def all_spectras(fpath='SpectraTa.csv', ofile='GammaPeaks.txt'):
    file = su.get_data(fpath)
    width = 1
    thr = 1000
    peak_lib = {}
    open(ofile, "w").close()
    f = open(ofile, "a")
    for i, specname in enumerate(file.spectras[1:]):

        calibration = float(file.slope[i])
        channels, counts = su.load_data(specname, file)
        energy = [i*calibration for i in channels]
        energy = np.asarray(energy)
        fig, axs = su.plotspec(energy, counts, title=specname,
                               xlabel='energy (keV)')

        # scipy find_peaks function
        peaks, props = find_peaks(counts, prominence=thr, width=width)
        k = []
        for i in peaks:
            try:
                guess = su.guess_fit(energy, counts, guess_0=energy[i],
                                     guess_2=1)
                xslice, yslice = su.slicer(energy, counts, guess)
                fit, popt = su.gauss_fit(xslice, yslice, guess=guess)
                area = simpson(fit, xslice)
                axs.plot(xslice, fit, 'r', ls='--')
                axs.annotate(f"{energy[i]:.2f} keV: {area:.2f} counts",
                             (energy[i], counts[i]))
                k.append(energy[i])
            except RuntimeError:
                print(specname, 'peak', energy[i])
        peak_lib.update({specname: k})
    for ele in peak_lib:
        f.write(f"{ele} {peak_lib[ele]}\n")
    f.close()


def iso_search(values):
    isos = dicts.iso_dict()
    for val in isos.values():
        peaks = []
        intensity = []
        for tuples in val:
            peaks.append(tuples[0])
            intensity.append(tuples[1])
            matcher(peaks, values)


def matcher(list1, list2, intensity):
    try:
        for i, ele in enumerate(list1):
            res = abs(ele - list2[i])
    except IndexError:
        pass
    return res
