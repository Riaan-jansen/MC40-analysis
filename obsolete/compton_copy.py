import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# import data from txt
filename = input('Please enter the name of the file : \t')
data = np.genfromtxt(filename, delimiter='\t', skip_header=3)

# extract x and y values from the data
x = data[:, 0]
ya = data[:, 1]

n = len(ya)
F = np.fft.fft(ya, n)
frequency = np.arange(n) / n    # freq. array
powers = F * np.conj(F) / n     # power series distribution
half_index = np.arange(1, np.floor(n/2), dtype=int)     
# floor of the scalar x is the largest integer i, such that i <= x.

plt.title('Power Spectrum Distribution. ' + filename + '.')
plt.xlabel('Frequency (KHz)')
plt.ylabel('Power (Amplitude)')
plt.plot(frequency[half_index], np.abs(powers[half_index]),
         label='Power Spectrum Distribution')
plt.legend()
plt.grid()
plt.show()

threshold = float(input('Input threshold power: '))
powers_index = powers > threshold   # array of 0 and 1 where good
powers_clean = powers * powers_index    # zeros all unnecessary powers
F_clean = powers_index * F

y_filtered = np.fft.ifft(F_clean)     # clean signal retrived

plt.title('Measurement count - Energy ' + filename + '(unfiltered)')
plt.plot(x, ya, 'r')
plt.ylabel('Counts')
plt.xlabel('Energy (Kev)')
plt.show()

ya = y_filtered

plt.plot(x, ya, 'steelblue')
plt.grid()
plt.title('Measurement count - Energy ' + filename)
plt.ylabel('Counts')
plt.xlabel('Energy (Kev)')
plt.show()

X = np.argsort(x)
x_new = np.zeros(len(X), dtype=float)

for i in range(0, len(X)):
    x_new[i] = x[X[i]]

plt.plot(X, ya)     # plots counts as function of indexes -
plt.grid()          # so that the data range can be selected for gaussian
plt.show()

lower_lim = input('Lower index limit for Gaussian fitting: ')
upper_lim = input('Upper index limit for Gaussian fitting: ')
x_gauss = x_new[int(lower_lim):int(upper_lim)]
y_gauss = ya[int(lower_lim):int(upper_lim)]


def func(x_gauss, *params):
    y = np.zeros(len(x_gauss), dtype=float)
    for i in range(0, len(params), 3):
        cntr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp(-((x_gauss - cntr) / wid)**2)
    return y


guess_0 = float(input('Input guess for Energy peak: '))
guess_1 = float(input('Input guess for Count max: '))
guess_2 = float(input('Input guess for FWHM: '))
guess = [guess_0, guess_1, guess_2]

popt, pcov = curve_fit(func, x_gauss, y_gauss, p0=guess, maxfev=8000)
fit = func(x_gauss, *popt)
sigma = np.abs(popt[2]) / 2.355
d_params = np.sqrt(np.diag(pcov))

print('Energy peak (central max) =\t %.3f ' % popt[0],
      '\nCount rate max =\t %.3f' % popt[1], '\nFWHM =\t %.3f' % popt[2])
print('Standard Deviation = \t', sigma)


plt.title('Compton Scatter . ' + filename + '. Count rate - Energy.')
plt.plot(x_new, ya, label='De-noised Data')
plt.plot(x_gauss, fit, 'r', ls='--', label='Gaussian Curve')
plt.ylabel('Counts')
plt.xlabel('Energy (Kev)')
plt.legend()
plt.grid()
plt.show()
