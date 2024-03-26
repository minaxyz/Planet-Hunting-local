import matplotlib.pyplot as plt
from data_analyser import DataAnalyser

periods = []
duration = []
for system in DataAnalyser():
    try:
        duration.append(system.getTransitDuration())
        periods.append(system.getOrbitalPeriod())
    except Exception:
        pass

plt.xlabel("Orbital Period")
plt.ylabel("Transit Duration")
plt.plot(periods, duration, '.')
plt.show()