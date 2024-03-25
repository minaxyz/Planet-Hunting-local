import matplotlib.pyplot as plt
from data_analyser import DataAnalyser

radii = []
for system in DataAnalyser():
    radii.append(system.getPlanetaryRadius())

plt.hist(radii)
plt.show()