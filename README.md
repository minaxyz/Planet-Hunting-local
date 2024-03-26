# Planet Hunting with Python

## Introduction

This is a repository that builds on provided code from QMUL. This is a module for analysing data from NASA's Kepler spacecraft.

## What this Repository contains
### data_analyser.DataAnalyser
Facilitates the extraction of system parameters from stellar system data.

It allows for the following system parameters to be extracted:
  - Stellar radius - The mean volumetric radius of the star -- Units: `Solar Radii`.
  - Stellar mass - The mass of the star -- Units: `Solar Masses`.
  - Orbital period - The time interval between transits -- Units: `Days`.
  - Phase - The time at which the first transit is detected -- Units: `Days`.
  - transit duration - The time the planet spends partially blocking the star's light -- Units: `Days`.
  - Planetary radius - The radius of the planet -- Units: `Earth Radii`.
  - Semi-major axis - The farthest distance between the planet and the star -- Units: `AU`.
  - Impact parameter - The perpendicular distance between the orbit the star's centre, expressed as a ratio of 
  the star's radius -- Units: `Solar radius ratio`.
  - Orbital inclination - The angle between the star's and the planet's orbiting plane -- Units: `degrees`.

### data_analyser.TransitDetector
Facilitates the detection of transits in flux against time data from stellar system.

### data_analyser.PhaseFoldedTransitModel
Creates a model for the phase-folded time-sorted transit data using polynomial interpolation.

### data_handler.AbstractDataHandler
Abstract data handler class for the creation of data handlers which are compatible with data analyser.

### data_handler.LocalDataHandler
Built-in fetching of data from specified local files which follow the old or new nasa exoplanet archive data format.

### formulas.py
A module containing all the formulas for the calculation of system parameters.

## Example Code
Plotting a histogram of planetary radii.
```
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
```

Plotting a graph of transit duration against orbital period.
```
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
```
