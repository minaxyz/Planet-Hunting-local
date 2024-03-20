from data_handler import LocalDataHandler
from data_analyser import DataAnalyser
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

G = 6.6743e-11 # The universal gravitational constant (m^3 kg^-1 s^-2).
SOLAR_MASS = 1.9885e30 # The volumetric mean radius of the sun (m).
SOLAR_RADIUS = 6.957e8 # The mass of the sun (kg).
EARTH_RADIUS = 6.378e6

''' Something Szymon did?? ''' 

# peakFlux = []
# stellar_radii = []
# for dataAnalyser in DataAnalyser():
#     peakFlux.append(dataAnalyser.getModel().getPeak())
#     stellar_radii.append(dataAnalyser.mass*6.957e+8)


# print(peakFlux)
# print(stellar_radii)
# planetary_radii = []
# for i in range(len(peakFlux)):
#     planetary_radii.append(round((stellar_radii[i]*(peakFlux[i]**0.5))/6.378e+6, 3))

# for i in planetary_radii:
#     planetary_radii[planetary_radii.index(i)] = round(i,2)

# print(planetary_radii)
# DataAnalyser('kplr002853093').plot('pm')
# print(DataAnalyser('kplr002853093').mass)
# print(DataAnalyser('kplr002853093').getModel().getPeak())

# def PlanetaryRadius(solar_radius, flux):
#     return (solar_radius*SOLAR_RADIUS*(abs(flux)**(1/2)))/EARTH_RADIUS

# print(PlanetaryRadius(1.2, 0.0045709))

''' End of Szymon's stuff '''

''' Histogram function '''

# with open('planet_radii.txt', 'r') as f:
#     data = [float(line.strip()) for line in f]
# hist, bins = np.histogram(data, bins='auto')
# plt.hist(bins[:-1], bins, weights=hist)
# plt.show()

''' End of histogram function '''