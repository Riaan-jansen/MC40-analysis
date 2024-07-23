'''file containing all the functions needed for spec_analysis notebook'''
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class file_data():
    '''class to contain data from csv file imported from Cambio.
    contains: spectras, rtimes, ltimes, channels, counts, and fdata'''
    def __init__(self):
        self.spectras = []
        self.rtimes = []
        self.ltimes = []
        self.channels = []
        self.counts = []
        self.fdata = None
        self.slope = None


class decay_data():
    def __init__(self):
        self.energy = []
        self.intensity = []


class functions():
    '''contains all the functions used for fitting and plotting. 
    contains gauss, exp growth and decay.'''
    def gauss(x, *params):
        '''gaussian function for curve fitting. parameters are:
        mean, height, FWHM. To normalise: A = 1/sqrt(2pi*wid**2).'''
        y = np.empty_like(x, dtype=float)
        for i in range(0, len(params), 3):
            cntr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            y = amp * np.exp(-((x - cntr) / wid)**2)
        return y

    def maxwell(x, *params):
        y = np.empty_like(x, dtype=float)

        for i in range(0, len(params), 3):
            A = params[i]
            B = params[i+1]
            C = params[i+2]
            y = A*(x**2)*np.exp(-((x - C) * B)**2)
        return y

    def expo(x, *params):
        '''given x (list/array), params = amplitude and lambda factor.
        returns array y.'''
        y = np.empty_like(x, dtype=float)
        for i in range(0, len(params), 3):
            A = params[i]
            B = params[i+1]
            C = params[i+2]
            y = A * x * np.exp(np.abs(x - C) * B)
        return y

    def series(x, *params):
        y = np.empty_like(x, dtype=float)
        for i in range(0, len(params), 3):
            a0 = params[i]
            a1 = params[i+1]
            a2 = params[i+2]
            y = a0 + a1*(1/x) + a2*(1/x)*2
        return y

    def yrs2sec(x):
        x = x * 365.25 * 24 * 60 * 60
        return x

    def surface_area(r):
        area = 4*np.pi*(r**2)
        return area


class dicts():
    '''class to store all dictionaries. g_dict is gamma peak energy
    and intensity. t_dict is half-life, activity(+error). cal_dict contains
    the channels for some spectra peaks. unit_dict contains mm and cm.'''
    def g_dict():
        '''dictionary containing the energy of gamma peaks of common
        calibration sources. e.g. dict = g_dict(); dict.keys() -> isotopes,
        values in tuple, [0] -> erg, [1] -> intensity %.'''
        gammas = {
            # values from lunds univ database
            '60Co': [(1173.2, 0.99974), (1332.5, 0.99986)],
            '133Ba': [(53.161, 0.022), (80.997, 0.3406), (276.40, 0.0716),
                      (302.85, 0.1833), (356.02, 0.6205), (383.85, 0.0894)],
            '152Eu': [(121.78, 0.2858), (244.70, 0.07583), (344.27, 0.2654),
                      (778.90, 0.1294), (964.08, 0.1461), (1085.9, 0.1021),
                      (1112.1, 0.1364), (1408.0, 0.2101)],
            '241Am': [(26.345, 0.024), (59.541, 0.3594)]
            # 'Background': [511, 583, 2611,   # Thallium-208 (Th)
            #                609, 1120, 1764,  # Bismuth-214 (U)
            #                911, 969,         # Actinium-228 (Th)
            #                239,              # Lead-212 (Th)
            #                1461]             # Potassium-40
        }
        return gammas

    def t_dict():
        '''dict containing the half-life of calibration sources in seconds,
        the activity at time of recording in Bq (1/s), and uncertainty.'''
        # yrs2sec = lambda x : x*365.25*24*60*60
        halflife = {
            '60Co': [functions.yrs2sec(5.2714), 39058, 1733],
            '133Ba': [functions.yrs2sec(10.51), 12981, 618],
            '152Eu': [functions.yrs2sec(13.537), 49095, 2455],
            '241Am': [functions.yrs2sec(432.2), 367961, 18398]
        }
        return halflife

    def cal_dict():
        '''contains calibrations specific to a file. instead of peak finder'''
        cals = {
            '23cm_60Co': (3440, 0),
            '23cm_152Eu': (4130, 7),
            '82mm_New_152Eu': (4090, 7)
        }
        return cals

    def unit_dict():
        units = {'cm': 1, 'mm': 0.1}
        return units

    def iso_dict(fpath='emissions_dat/'):
        W181 = read_emissions(fpath+'W181.csv')
        W180 = read_emissions(fpath+'W180.csv')
        W179 = read_emissions(fpath+'W179.csv')
        Hf178m1 = read_emissions(fpath+'Hf178m1.csv')
        Hf178m2 = read_emissions(fpath+'Hf178m2.csv')
        isos = {
            '181W': W181,
            '180W': W180,
            '179W': W179,
            '178Hfm1': Hf178m1,
            '178Hfm2': Hf178m2
        }
        return isos


def read_emissions(fpath):
    tuples = []
    with open(fpath) as f:
        lines = f.read().splitlines()
        lines = [row.split(',') for row in lines[1:]]
        for row in lines:
            if row[4] == 'G':
                tuples.append((float(row[0]), float(row[2])))
    return tuples


def read_emissions_lnhb(fpath):
    '''obsolete-
    given filepath, reads .txt and returns energy, intensity tuple.
    for LNHB lab database .lara emissions files only.'''
    with open(fpath) as f:
        lines = f.read()
        lines = lines.replace(' ', '')
        lines = lines.splitlines()
        lines = lines[14:-1]
        lines = [row.split(';') for row in lines]
        erg = [float(row[0]) for row in lines]
        int = [float(row[2]) for row in lines]
        if len(erg) != len(int):
            raise ValueError('lists not the same length in emissions file')
        tuple = [(ele, int[i]) for i, ele in enumerate(erg)]
    return tuple


def get_data(filename):
    '''reads data from Cambio into list.
    returns file object.'''
    file = file_data()

    file.fdata = reader(filename)
    file.spectras = file.fdata[0]

    file.rtimes = file.fdata[4]
    file.rtimes = file.rtimes[1:]
    file.ltimes = file.fdata[5]
    file.ltimes = file.ltimes[1:]
    file.rtimes = [float(i) for i in file.rtimes]
    file.ltimes = [float(i) for i in file.ltimes]
    file.slope = file.fdata[7]
    file.slope = file.slope[1:]

    file.channels = []
    file.counts = []
    return file


def load_data(specname, file, headers=13):
    '''given spectra name and read file, extracts channel and count.'''
    for i, name in enumerate(file.spectras):
        if specname == name:
            index = i
    file.channels = []
    file.counts = []
    for i, line in enumerate(file.fdata):
        if i >= headers:
            file.channels.append(float(line[0]))
            file.counts.append(float(line[index]))
    return file.channels, file.counts


def gauss_fit(x, y, guess=None):
    '''fits a gaussian given x, y and and guess.
    outputs: fit (y values), and optimised parameters.'''
    popt, cov = curve_fit(functions.gauss, x, y, p0=guess, maxfev=8000)
    fit = functions.gauss(x, *popt)
    return fit, popt


def guess_fit(x, y, guess_0=None, guess_2=7):
    '''guess x value for y peak, and width of peak (Default=7) to minimise
    curve_fit run time, and point to relevant peak. returns guess:
    list = [mean, amplitude, width].'''
    if guess_0 is None:
        guess_0 = float(input('Input guess for channel peak: '))
        guess_index = x.index(guess_0)
    else:
        guess_array = np.empty_like(x)
        for i, ele in enumerate(x):
            guess_array[i] = np.abs(ele - guess_0)
            guess_index = guess_array.argmin()

    guess_1 = y[guess_index]
    guess = [guess_0, guess_1, guess_2]
    return guess


def slicer(x, y, guess):
    '''takes x, y and guess for center of gaussian, and slices the xdata
    around it. returns x and y data where anywhere outside center
    +/- range = 0.'''
    def chopper(a, a_slice):
        '''obsolete -
        if elements of a are not in a_slice they are = 0. returns a'''
        for i, ele in enumerate(a):
            if ele not in a_slice:
                a[i] = 0
        return a
    width = int(guess[2])-2
    array = np.empty_like(x)
    for i, ele in enumerate(x):
        array[i] = np.abs(ele - guess[0])
        close_index = array.argmin()
    x_slice = x[close_index-width:close_index+width]
    y_slice = y[close_index-width:close_index+width]
    return x_slice, y_slice


def fft_cleaning(xdata, ydata):
    '''function to denoise data if needed for gaussian fitting'''
    n = len(ydata)
    F = np.fft.fft(ydata, n)
    # freq. array
    frequency = np.arange(n) / n
    # power series distribution
    powers = F * np.conj(F) / n
    # floor of the scalar x is the largest integer i such that i <= x.
    half_array = np.arange(1, np.floor(n/2), dtype=int)

    plt.title('Power Spectrum Distribution.')
    plt.xlabel('Frequency (KHz)')
    plt.ylabel('Power (Amplitude)')
    plt.plot(frequency[half_array], np.abs(powers[half_array]),
             label='Power Spectrum Distribution')
    plt.legend()
    plt.grid(which='both')
    plt.minorticks_on()
    plt.show(True)

    threshold_value = 0.001
    threshold = threshold_value * np.max(np.abs(powers))
    # array of 0 and 1 for below/above threshold values
    powers_index = powers > threshold
    # powers_clean = powers * powers_index
    # zeros all unnecessary powers
    F_clean = powers_index * F
    y_filtered = np.fft.ifft(F_clean)

    plt.plot(xdata, y_filtered, 'steelblue')
    plt.grid(which='both')
    plt.minorticks_on()
    plt.title('Filtered Data.')
    plt.ylabel('y-axis')
    plt.xlabel('x-axis')
    plt.show()

    return y_filtered


def surface_ratio(specname, r=3):
    '''calculates the ratio of detector surface to emitting area (4 pi r^2)
    uses spectrum filename, radius of detector (default=3). returns ratio'''
    units = dicts.unit_dict()
    for ele in units.keys():
        if ele in specname:
            unit_factor = units[ele]
            specname = specname.split(ele)
            R = float(specname[0])
        else:
            pass
    R = R * unit_factor
    area_det = functions.surface_area(r)
    area_emit = functions.surface_area(R)
    ratio = float(area_det/area_emit)
    return ratio


def plotspec(x, y, title=None, xlabel='x-axis'):
    '''basic plotting function. args = x, y; kwargs = title, xlabel'''
    fig, axs = plt.subplots(figsize=(10, 6))
    fig.suptitle(title + ' Spectrum.')
    fig.supylabel('Counts')
    fig.supxlabel(xlabel)
    axs.grid(which='major', lw=1.3)
    axs.grid(which='minor', linewidth=0.5)
    axs.minorticks_on()
    axs.plot(x, y, color='blue', label='Spectrum')
    fig.legend()
    return fig, axs


def reader(filename):
    '''reads data from csv format file into a list of lists'''
    with open(filename) as f:
        lines = f.read().splitlines()
        lines = [row.split(',') for row in lines]
    return lines


def efficiency_calc(N_counts: float, time: float, Yield: float,
                    isotope: str, specname: str, absolute=True):
    '''calculates the detector efficiency given total counts, real time,
    initial activity and source half-life. returns efficiency for
    the given energy.'''
    db = dicts.t_dict()
    for ele in db.keys():
        if ele == isotope:
            t_2, A_0, A_0err = db[ele]
    activity = A_0 * np.exp(-np.log(2) * time / t_2)

    N_emit = activity * Yield * time

    area_ratio = surface_ratio(specname)
    efficiency = float(N_counts / (N_emit * area_ratio))

    if absolute is True:
        efficiency = float(N_counts / N_emit)

    return efficiency*100


def plot_eff(energies, efficiencies):
    '''function to plot efficiency - energy'''
    fig, axs = plt.subplots(figsize=(10, 6))
    fig.suptitle('Efficiency - Energy plot for semiconductor detector.')
    fig.supylabel('Efficiency (%)')
    fig.supxlabel('Energy (keV)')
    axs.grid(which='major', lw=1.3)
    axs.grid(which='minor', linewidth=0.5)
    axs.minorticks_on()
    axs.scatter(energies, efficiencies, color='blue')
    fig.legend()
    return fig, axs


def unpacker(list):
    '''unpacks list of tuples and sorts in ascending order in x
    (first element of tuple). returns x and y (lists).'''
    x = []
    yi = []
    y = []
    for tuple in list:
        x.append(tuple[0])
        yi.append(tuple[1])
    sortedx = np.argsort(x)
    x.sort()
    for i in sortedx:
        y.append(yi[int(i)])
    return x, y


def eff_curve(energy, efficiency, d='long'):
    if d == 'long':
        guess = [0.2, -0.00568, 53]  # guesses come from evaluating
    else:                            # data at points to find exp
        guess = [1.2, -0.00568, 53]  # coefficients
    eff_array = np.array(efficiency, dtype=float)
    index = (eff_array.argmax()) + 1
    erg1 = energy[:index]
    eff1 = efficiency[:index]
    erg2 = energy[index:]
    eff2 = efficiency[index:]

    popt2, cov2 = curve_fit(functions.expo, erg2, eff2, p0=guess)
    fit2 = functions.expo(erg2, *popt2)
    fig, ax = plot_eff(energy, efficiency)
    if d == 'long':
        c = 'red'
    else:
        c = 'orange'

    ax.plot(erg1, eff1, color=c, label='efficiency curve (' + d + ')')
    ax.plot(erg2, fit2, color=c)
    fig.legend()
    return fig, ax
