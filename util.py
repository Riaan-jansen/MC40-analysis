import spec_utilities as su
import numpy as np
import matplotlib.pyplot as plt

def surface_area(specname, r=3):
    '''calculates the ratio of detector surface to emitting area (4 pi r^2)
    uses spectrum filename, radius of detector (default=3). returns ratio'''
    units = {'cm': 1, 'mm': 0.1}
    for ele in units.keys():
        if ele in specname:
            unit_factor = units[ele]
            specname = specname.split(ele)
            R = float(specname[0])
        else:
            pass
    R = R * unit_factor
    area_det = r**2
    area_emit = 4*(R**2)
    ratio = float(area_det/area_emit)
    return ratio


# filename = 'Spectra.csv'
# file = su.get_data(filename)

# for j, specname in enumerate(file.spectras[1:-1]):

#     tuple = surface_area(specname)
#     print(tuple)

x = np.arange(-np.pi, np.pi, np.pi/100)
y = np.cos(x)

plt.plot(x, y)
plt.show()
