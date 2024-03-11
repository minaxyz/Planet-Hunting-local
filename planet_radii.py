from data_handler import LocalDataHandler
from data_analyser import DataAnalyser
from time import perf_counter
import matplotlib.pyplot as plt
import numpy as np

# peakFlux = []
# stellar_radii = []
# for dataAnalyser in DataAnalyser():
#     peakFlux.append(dataAnalyser.getModel().getPeak())
#     stellar_radii.append(dataAnalyser.mass*6.957e+8)


# # print(peakFlux)
# # print(stellar_radii)
# planetary_radii = []
# for i in range(len(peakFlux)):
#     planetary_radii.append(round((stellar_radii[i]*(peakFlux[i]**0.5))/6.378e+6, 3))

# for i in planetary_radii:
#     planetary_radii[planetary_radii.index(i)] = round(i,2)

# print(planetary_radii)

print(DataAnalyser('kplr002853093').getModel().getPeak())
print(LocalDataHandler())

def getPlanetaryRadius(self):
        r = self.radius*6.957e+8
        flux = DataAnalyser().getModel().getPeak()
        pr = r*(flux**(1/2))/6.378e+6
        return pr