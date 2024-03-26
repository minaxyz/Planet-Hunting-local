import matplotlib.pyplot as plt
from data_analyser import DataAnalyser

radii = []
for system in DataAnalyser():
    try:
        radii.append(system.getPlanetaryRadius())
    except Exception:
        pass

plt.xlabel("Planetary Radius")
plt.ylabel("Frequency")
plt.hist(radii)
plt.show()