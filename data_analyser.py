import matplotlib.pyplot as plt
import formulas
import numpy as np
import inspect
from numpy.polynomial.polynomial import Polynomial, polyval, polyder, polyroots
from math import floor
import inspect

from data_handler import AbstractDataHandler, LocalDataHandler
from utils import estimatePeriodicSignal

##Transit Detection Constants
FLUX_SAMPLES = 1000
SIGNIFICANCE_LEVEL = 0.05 #Determines the transit threshold
MINIMUM_PERIOD = 0.241842217 #Minimum period used to determine the transit detector uniform convolution factor. Sourced from the NASA exoplanet archive.
#Anomalous Region Constants
ANOMALOUS_REGION_SEARCH_OFFSET = 1 #The area to search for anomalous points (in MINIMUM_PERIODS).
ANOMALOUS_CONCENTRATION_THRESHOLD = 0.1 #The concentration of anomalous points for it to be considered.
TIME_TO_SEARCH_FOR_ANOMALOUS_REGIONS = 10

##Transit Analysis Constants
#Transit Calibration Constants
TRANSIT_SAMPLES = 30 #Number of samples to complete per iteration.
CALIBRATION_SEARCH_OFFSET = 0.1 #Determines the interval to search for the transit.
TRANSIT_THRESHOLD_ITERATION_SCALING = 0.75 #The scaling applied to the threshold after each iteration. Must be between 0 and 1.
MINIMUM_TRANSIT_THRESHOLD = 0.5 #Determines the maximum number of recursion of calibration (max recursion depth = ceil(log_TRANSIT_THRESHOLD_ITERATION_SCALING(MINIMUM_TRANSIT_THRESHOLD))).
"""Determines the strictness in validating a transit (Higher -> Higher strictness), must be between 0 and 1, higher than ACCEPTANCE_LEVEL.
Currently, accepts a transit after 3 passes and 0 fails, or accepted if above 80% passes. (#Passes no fails to accept = (1/1 - ACCEPTANCE_LEVEL) - 2))"""
ACCEPTANCE_LEVEL = 0.8
"""Determines the strictness in rejecting a transit (Lower -> Higher strictness), must be between 0 and 1, lower than REJECTION_LEVEL.
Currently, rejects a transit after 0 passes and 3 fails, or rejected if below 20% passes. (#Fails no passes to reject = (1/REJECTION_LEVEL) - 2"""
REJECTION_LEVEL = 0.2
TRANSIT_CHARACTERISTICS_CLOSENESS_LEVEL = 0.025 #Determines if a period and phase of a transit pair partains to the current period and phase found.
##Orbital Period Calculation Constants
PERIOD_CALC_SEARCH_OFFSET = 0.2 #Determines the interval to search for the transit (Fractions of the period).
RECALIBRATION_FACTOR = 2 #Determines how often to recalibrate (Recalibration occurs every log_RECALIBRATION_FACTOR(NTransits) times).

def plot(*plots):
    """Uses matplotlib to plot a time aganinst flux graph.

    ----------
    Arguments:
        plot (times, flux):
            times (np.array) - Contains the times for the corresponding flux recordings.

            flux (np.array) - Contains the recordings of the intensity of the star's light.
        or
        *plots (arraylike):
            *plots - An array of individual plots (`(plot1, plot2,.. plotN), where plotK = (timesK, fluxK)`).
    """
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("Flux")
    if isinstance(plots[0][0], float):
        plt.plot(*plots, '.', markersize=2)
    else:
        for data in plots:
            plt.plot(*data, '.', markersize=2)
    plt.show()

def histogram(flux):
    """Uses matplotlib to plot a histogram of the flux values. The histogram contains 100 buckets.
    The y-axis is the percentage of the flux values which an individual bucket occupies.

    ----------
    Arguments:
        flux (np.array) - The array containing the flux values to form a histogram from.
    """
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("Proportion")
    plt.hist(flux, 100, weights=np.full(len(flux),1/len(flux),dtype=float))
    plt.show()

def phaseFold(times, flux, period, phase):
    """Returns an array of phase-folded times and flux, sorted by phase-folded times.

    ----------
    Arguments:
        times (np.array) - Contains the times for the corresponding flux recordings.

        flux (np.array) - Contains the recordings of the intensity of a star's light intensity over time.

        period (float) - The period of the transits.

        phase (float) - The time at which the peak of any of the transits occurs.
    
    ----------
    Returns:
        times (np.array) - Contains the times for the corresponding flux recordings.

        flux (np.array) - Contains the recordings of the intensity of a star's light intensity over time.
    """
    phaseFoldedTimes = ((times - phase + period/2) % period) - period/2
    sort = np.argsort(phaseFoldedTimes)
    return phaseFoldedTimes[sort], flux[sort]

class TransitDetector():
    def __init__(self, times, flux, searchMode=True):
        """
        Arguments:
            times (np.array) - Contains the times for the corresponding flux recordings.

            flux (np.array) - Contains the recordings of the intensity of a star's light intensity over time.

            searchMode (bool) - When True, reduces noise to enhance transits, making them easier to detect. 
            If quality needs to be preserved, it is recommended to set this to False (True by default).
        """
        self.searchMode = searchMode
        self.times = times
        self.size = len(self.times)
        self.dt = (self.times[-1] - self.times[0])/self.size
        #Applies a convolution to the flux which reduces noise and accentuates transits.
        self.uniformConvolutionFactor = floor((0.5 if searchMode else 0.05)*MINIMUM_PERIOD/self.dt) or 1
        
        #Removal of outliers if searchMode is off.
        if not searchMode: 
            self.convolvedFlux = flux
            self.medianFlux = None
            self.transitThreshold = None
            anomalousThreshold = -self.getTransitThreshold()
            anomalousIndexes = [i for i, flux in enumerate(flux) if flux > anomalousThreshold]
            self.times = np.delete(self.times, anomalousIndexes)
            flux = np.delete(flux, anomalousIndexes)
        
        self.convolvedFlux = np.convolve(flux, np.full(self.uniformConvolutionFactor, 1/self.uniformConvolutionFactor, dtype=float),'same')
        
        self.anomalousRegions = [] #Even indexes represent the start of an anomalous region, and its next odd index represents the end of the anomaly.
        self.medianFlux = None #The median flux, calculated from the first FLUX_SAMPLES flux values.
        self.transitThreshold = None #The threshold below which transits should occur. Implementation details in self.getTransitThreshold.
        
        self.start = self.times[0]
        self.end = self.times[-1]

    def __findTime(self, time):
        return min(self.times.searchsorted(time), self.size - 1)

    def findTransit(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None):
        """Returns the time at which the first transit was detected, None if no transit was detected in the given time bounds.

        ----------
        Arguments:
            timeStart (float) - The time at which to start searching for transits (The first value of the time array by default, last if reverse).
            
            timeEnd (float) - The time at which to stop searching for transits (The last value of the time array by default, first if reverse).

            reverse (bool) - The direction in which to search for transits (False by default).

            transitThreshold (float) - The flux level below which a transit is detected, relative to the standard transit threshold 
            (see `getTransitThreshold()` for more information) (1 by default).

        ----------
        Returns:
            time (float|None) - The time at which the transit was detected, or None if no transit was detected.
        """
        self.getTransitThreshold()
        i = self.__findTransit(
            self.__findTime(timeStart) if timeStart is not None else (self.size if reverse else 0),
            self.__findTime(timeEnd) if timeEnd is not None else (-1 if reverse else self.size), 
            reverse, 
            self.transitThreshold if transitThreshold is None else self.medianFlux + transitThreshold*(self.transitThreshold - self.medianFlux))
        return i and self.times[i]

    def findTransitPeak(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None, searchMode=True):
        """Returns the time of the peak of the first transit detected, None if no transit was detected in the given time bounds.

        ----------
        Arguments:
            timeStart (float) - The time at which to start searching for transits (The first value of the time array by default, last if reverse).
            
            timeEnd (float) - The time at which to stop searching for transits (The last value of the time array by default, first if reverse).

            reverse (bool) - The direction in which to search for transits (False by default).

            transitThreshold (float) - The flux level below which a transit is detected, relative to the standard transit threshold 
            (see `getTransitThreshold()` for more information) (1 by default).

            searchMode (bool) - If true, will search the given time bounds sequantially from timestart to timeEnd, else will search and return the 
            transit closest to the midpoint of timeStart and timeEnd.

        ----------
        Returns:
            time (float|None) - The time of the peak of the first transit detected, or None if no transit was detected.
        """
        return self.findTransitBounds(timeStart, timeEnd, reverse, transitThreshold, searchMode)[1]
    
    def findTransitBounds(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None, searchMode=True):
        """Returns the start, peak, and end time of the peak of the first transit detected, 
        None if no transit was detected in the given time bounds.

        ----------
        Arguments:
            timeStart (float) - The time at which to start searching for transits (The first value of the time array by default, last if reverse).
            
            timeEnd (float) - The time at which to stop searching for transits (The last value of the time array by default, first if reverse).

            reverse (bool) - The direction in which to search for transits (False by default).

            transitThreshold (float) - The flux level below which a transit is detected, relative to the standard transit threshold 
            (see `getTransitThreshold()` for more information) (1 by default).

            searchMode (bool) - If true, will search the given time bounds sequantially from timestart to timeEnd, else will search and return the 
            transit closest to the midpoint of timeStart and timeEnd.

        ----------
        Returns:
            tuple (start, peak, end):
                start (float) - The time at which the transit is detected.
                
                peak (float) - The time at which the transit reaches its maximum flux.
                
                end (float) - The time at which the transit ends.
            or (None, None, None) if no transit was detected in the given time bounds.
        """
        self.getTransitThreshold()
        if searchMode:
            bounds = self.__findTransitBounds(
                self.__findTime(timeStart) if timeStart is not None else (self.size if reverse else 0),
                self.__findTime(timeEnd) if timeEnd is not None else (-1 if reverse else self.size), 
                reverse, 
                self.transitThreshold if transitThreshold is None else self.medianFlux + transitThreshold*(self.transitThreshold - self.medianFlux))    
        else:
            midPoint = ((timeStart + timeEnd)/2 if timeEnd is not None else timeStart) if timeStart is not None else (self.start + self.end)/2
            searchOffset = (timeEnd - timeStart)/2 if timeStart is not None and timeEnd is not None else (self.end - self.start)
            leftTransitBounds = self.findTransitBounds(midPoint, midPoint + searchOffset, False, transitThreshold)
            rightTransitBounds = self.findTransitBounds(midPoint, midPoint - searchOffset, True, transitThreshold)
            if leftTransitBounds[1] is None:
                return rightTransitBounds
            elif rightTransitBounds[1] is None:
                return leftTransitBounds
            else:
                return leftTransitBounds if leftTransitBounds[0] - midPoint < midPoint - rightTransitBounds[2]  else rightTransitBounds
        return (None, None, None) if bounds is None else tuple(self.times[x] for x in bounds)

    def __findTransitBounds(self, start, end, reverse, transitThreshold):
        """
        Returns the time index of a transit's start, peak and end.

        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            revese (bool) - Reverses the direction of the search if True, by default the timestep is the same as the time array.

        ----------
        Returns:
            tuple (start, peak, end)|None:
                start (int) - The time index at which the transit is detected.
                
                peak (int) - The time index at which the transit reaches its maximum flux.
                
                end (int) - The time index at which the transit ends.
            or None if no transit was detected in the given time bounds.
        """
        START = start
        END = end
        #Ensuring that the transit is not anomalous.
        while (start := self.__findTransit(start, end, reverse, transitThreshold)) != (nregion := self.__findNormalRegion(start, reverse)) and self.searchMode:
            start = nregion
        #Returning if None i.e. transit is not found in the given bound.
        if start is None or (start < end if reverse else start > end):
            return None
        #Decrement start, until the start of the transit is found.
        for start in range(start, -1, -1):
            if self.convolvedFlux[start] > transitThreshold:
                break
        #Iterate forwards, keeping count of the maximum flux and its location, until the end of the transit is reached.
        fluxLow = self.convolvedFlux[start]
        iLow = start
        i = start
        for i in range(start + 1, self.size, 1):
            if (transitBound := self.convolvedFlux[i]) < fluxLow:
                fluxLow, iLow = transitBound, i
            elif self.convolvedFlux[i] > transitThreshold:
                return (start, iLow, i) or None if (END <= iLow <= START if reverse else START <= iLow <= END) else None
        return start, iLow, i

    def __findTransit(self, start, end, reverse, transitThreshold):
        """
        Searches for the start of a transit. For a transit to be detected, the convolved flux value must be below the transit bound.
        
        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            revese (bool) - Reverses the direction of the search if True, by default the timestep is the same as the time array.

        ----------
        Returns:
            transit bound (int|None) - The time index at which the transit is detected. None is returned if a transit is not found.
        """
        if start is None:
            return None
        start = min(max(start, 0), self.size - 1)
        for i in range(start, end, -1 if reverse else 1):
            if self.convolvedFlux[i] < transitThreshold:
                return i
        return None
    
    def getTransitThreshold(self):
        """
        Returns an approximate flux value below which transits should occur.

        ----------
        Returns:
            flux value (float) -- Defined as `Median convolved flux value - 3(SIGNIFICANCE_LEVEL from the highest convolved flux value)`
            from a sample of the first SAMPLE values.
        """
        if self.transitThreshold is None:
            sortedSamples = sorted(self.convolvedFlux[:FLUX_SAMPLES])
            self.medianFlux = sortedSamples[floor(0.5*FLUX_SAMPLES)]
            self.transitThreshold = self.medianFlux + 3*(self.medianFlux - sortedSamples[floor((1 - SIGNIFICANCE_LEVEL)*FLUX_SAMPLES)])
        return self.transitThreshold

    def __findNormalRegion(self, start, reverse):
        if start is None:
            return None
        #If anomalous region already discovered, return the end bound of the region.
        
        if ((anomalousRegion := np.searchsorted(self.anomalousRegions, start)) & 1) == 1:
            return self.anomalousRegions[anomalousRegion - 1] if reverse else self.anomalousRegions[anomalousRegion]
        #Checking if the region within searchOffset of the index is anomalous.
        direction = -1 if reverse else 1
        timeSearchOffset = direction*MINIMUM_PERIOD*ANOMALOUS_REGION_SEARCH_OFFSET
        anomalyThreshold = -self.transitThreshold
        anomalyCount = 0
        firstAnomaly = None
        for i in range(self.__findTime(self.times[start] - timeSearchOffset), self.__findTime(self.times[start] + timeSearchOffset), direction):
            if self.convolvedFlux[i] > anomalyThreshold:
                if anomalyCount == 0:
                    firstAnomaly = i - direction
                elif anomalyCount >= 10:
                    break
                anomalyCount += 1
        #If the number of anomalies is more than 3, the region is anomalous.
        if firstAnomaly is None or anomalyCount < ANOMALOUS_CONCENTRATION_THRESHOLD*2*timeSearchOffset/self.dt:
            return start
        i = start
        step = direction*floor(TIME_TO_SEARCH_FOR_ANOMALOUS_REGIONS/self.dt)
        end = i + step
        while i != end and 0 < end < self.size:
            if self.convolvedFlux[i] > anomalyThreshold:
                end = i + step
            i += direction
        anomalyEnd = end - step
        self.__addAnomalousRegion(*((anomalyEnd, firstAnomaly) if reverse else (firstAnomaly, anomalyEnd)))
        return anomalyEnd
    
    def addAnomalyRegion(self, timeStart, timeEnd):
        """Marks a region as anomalous. Transits are ignored within anomalous regions.

        ----------
        Arguments:
            start (float) - The time of the start of the anomalous region.
            end (float) - The time to end of the anomalous region.
        """
        self.__addAnomalyRegion(self.__findTime(timeStart), self.findTime(timeEnd))

    def __addAnomalousRegion(self, start, end):
        searchOffset = floor(MINIMUM_PERIOD/(2*self.dt)) + 1
        start -= searchOffset
        end += searchOffset
        startI = np.searchsorted(self.anomalousRegions, start)
        endI = np.searchsorted(self.anomalousRegions, end)
        if startI == len(self.anomalousRegions) or endI == 0: #Needs to be placed at the start or end of transit data.
            self.anomalousRegions[startI:startI] += [start, end]
        else:
            #Resolution of potentially overlapping regions.
            isBoundOutsideRegion = lambda boundI : (boundI & 1) == 0 #If bound in of an even index, it is outside a region.
            self.anomalousRegions[startI:endI] = [bound for bound, boundI in zip([start, end], [startI, endI]) if isBoundOutsideRegion(boundI)]
    
    def getAnomalousRegions(self):
        """Returns the anomalous regions (regions within which transits are ignored).

        ----------
        Returns:
            *anomalous regions (array[(float, float)]) - An array of the anomalous regions, where anomalous regions are represented by
            a tuple of the time at which they start and end.
        """
        return [[self.times[x] for x in self.anomalousRegions[i:i+2]] for i in range(0,len(self.anomalousRegions),2)]

    def getData(self):
        """Returns the anomalous regions.

        ----------
        Returns:
            tuple (times, convolved flux)
                times (np.array) - Contains the times for the corresponding flux recordings.

                convolved flux (np.array) - The flux recordings convolved to reduce noise, and so enhance transits.
        """
        return self.times, self.convolvedFlux

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = self.__findTime(key.start)
            stop = self.__findTime(key.stop)
            #TODO: Prevent values above the anomaly threshold from being returned (Interferes with the transit model).
            return self.times[start:stop], self.convolvedFlux[start:stop]
        else:
            i = self.__findTime(key)
            return self.times[i], self.convolvedFlux[i]

class DataAnalyser(AbstractDataHandler):
    """Class for the analysis of stellar systems.
    Allows for system parameters to be extrancted from the required data provided by an instance of a subclass of AbstractDataHandler.

    System parameters supported:
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
    """

    def __init__(self, dataID=None, highAccuracy=True, dataHandler:AbstractDataHandler=LocalDataHandler):
        """
        Arguments:
            dataID (str) - The ID of the data to load.

            highAccuracy (bool) - If True, calculates to the highest degree of accuracy at the cost of performance. True by default.
            
            dataHandler (AbstractDataHandler) - Class or instance to be used to fetch required data from 
            (times against flux arrays, stellar radius, stellar mass), must implement AbstractDataHandler methods. 
            The LocalDataHandler class is used by default.
        """
        self.dataHandler = dataHandler(dataID) if dataID else dataHandler if not inspect.isclass(dataHandler) else None
        self.dataID = dataID or self.dataHandler and self.dataHandler.dataID
        #Flux against Time data
        self.times, self.flux = self.dataHandler.getData() if self.dataHandler is not None else (None, None)


        self.phaseFoldedTimes, self.phaseFoldedFlux = None, None
        self.transits = TransitDetector(self.times, self.flux) if self.dataHandler is not None else None
        self.model = None
        #Stellar radius and mass
        self.radius, self.mass = (self.dataHandler.getRadius(), self.dataHandler.getMass()) if self.dataHandler is not None else (None, None)

        self.size = len(self.times) if self.dataHandler is not None else None
        self.period = None
        self.transitThreshold = None
        self.transitDuration = None
        self.phase = None

        self.HIGH_ACCURACY_FLAG = highAccuracy
        self.CALIBRATION_FLAG = False
        self.SAMPLING_FLAG = False

    def plot(self, plotType=""):
        """Uses matplotlib.pyplot to plot a graph using the data provided to DataAnalyser via the TransitModel.

        ----------
        Arguments:
            plotType (str) - Specifies the type of plot to plot. Must be one of the following options:
                Standard ("phase" | "phase folded" | "p") - Plots the time against flux data.

                Phase Folded ("phase" | "phase folded" | "p") - Plots the time, phase-folded by the orbital period calculated, against flux data.

                Convolted ("convolved" | "convolution" | "c") - Plots the time against the convolted flux data.

                Model ("model" | "m") - Plots the model's time against flux sampled at the data's time recording interval.

                Phase model ("phase model" | "pm") - Plots the both the model's and phase folded time against flux data.

                Convolved phase ("convolved phase" | "cp") - Plots the phase-folded time against the convolted flux data.

                Convolved phase model ("convolved phase model" | "cpm")- Plots both model and the phase-folded time against the convolted flux data.

                Histogram ("histogram" | "h") - Forms a histogram of fluxes in the flux data.

                Convolved histogram ("convolved histogram" | "ch") - Forms a histogram of fluxes in the convolved flux data.
        """
        match plotType.lower():
            case "normal" | "standard" | "n" | "s" | "":
                plot(self.getData())
            case "phase" | "phase folded" | "p":
                plot(self.getPhaseFoldedData())
            case "convolved" | "convolution" | "c":
                plot((self.times, self.transits.convolvedFlux))
            case "model" | "m":
                plot(self.getModel().getData())
            case "phase model" | "pm":
                plot(self.getPhaseFoldedData(), self.getModel().getData())
            case "convolved phase" | "cp":
                plot(self.getModel().transitDetector.getData())
            case "convolved phase model" | "cpm":
                plot(self.getModel().transitDetector.getData(), self.getModel().getData())
            case "histogram" | "h":
                histogram(self.flux)
            case "convolved histogram" | "ch":
                histogram(self.transits.convolvedFlux)
            case _:
                raise Exception("Invalid Plot Type: Plot not recognised.\nPlot Type Options: 'standard', 'phase folded', 'model', 'phase model'.")
        plt.show()

    def getPhaseFoldedData(self):
        """Returns the flux against time data of the stellar system.
        
        ----------
        Returns:
            tuple (phase-folded times, flux):
                phase-folded times (np.array) - Contains the phase-folded times for the corresponding flux recordings.
                Generated by `((times - phase + period/2) % period) - period/2`.

                flux (np.array) - Contains the recordings of the intensity of the star's light. Data must be normalised
                (i.e. flux must be a fraction of the mean).

            *Note that the arrays returned must not contain NaN values.
        """
        if self.phaseFoldedTimes is None:
            self.getOrbitalPeriod()
            self.phaseFoldedTimes, self.phaseFoldedFlux = phaseFold(self.times, self.flux, self.period, self.phase)
        return self.phaseFoldedTimes, self.phaseFoldedFlux
    
    def getModel(self):
        """Returns the instance of the model of the phase-folded transit.

        Returns:
            model (PhaseFoldedTransitModel) - The instace of the model of the phase-folded transit.
        """
        if self.model is None:
            self.model = PhaseFoldedTransitModel(*self.getPhaseFoldedData())
            self.transitDuration = self.model.max - self.model.min
        return self.model
    
    def getOrbitalPeriod(self):
        """Returns the lowest planetary orbital period. 

        Returns:
            orbital period (float) - The time the planet takes to complete and orbit around its star -- Units: `Days`.
        """
        if self.HIGH_ACCURACY_FLAG and not self.SAMPLING_FLAG:
            self.__sampleTransits()
        elif not self.CALIBRATION_FLAG:
            self.__calibrate()
        return self.period

    def getTransitDuration(self):
        """Returns the transit duration of the planet with the lowest orbital period. 

        Returns:
            transit duration (float) - The time the planet spends partially blocking the star's light -- Units: `Days`.
        """
        if self.HIGH_ACCURACY_FLAG and not self.model:
            self.getModel()
        if self.transitDuration is None:
            start, peak, end = self.transits.findTransitBounds(self.getPhase(), transitThreshold=self.transitThreshold)
            if start and end:
                self.transitDuration = max(end - start - 2*self.transits.dt*self.transits.uniformConvolutionFactor, self.transits.dt)
        return self.transitDuration

    def getPhase(self):
        """Returns the transit duration of the planet with the lowest orbital period. 

        Returns:
            transit duration (float) - The time the planet spends partially blocking the star's light -- Units: `Days`.
        """
        if self.phase is None:
            self.__calibrate()
        return self.phase
    
    def getPlanetaryRadius(self):
        """Returns the planetary radius of the planet with the lowest orbital period. 

        Returns:
            planetary radius (float) - The radius of the planet -- Units: `Earth Radii`.
        """
        return self.radius and self.getModel().getPeak() and formulas.planetaryRadius(self.radius, self.getModel().getPeak())
    
    def getSemiMajorAxis(self):
        """Returns the semi-major axis of the planet with the lowest orbital period. 

        Returns:
            semi-major axis (float) - The farthest distance between the planet's and the star's centres -- Units: `AU`.
        """
        if self.mass and self.getOrbitalPeriod():
            return self.mass and formulas.semiMajorAxis(self.mass, self.period)
        return None

    def getImpactParameter(self):
        """Returns the impact parameter of the planet with the lowest orbital period. 

        Returns:
            impact parameter (float) - The perpendicular distance between the orbit the star's centre, expressed as a ratio of the 
            star's radius -- Units: `Solar radius ratio`.
        """
        if self.radius and self.mass and self.getOrbitalPeriod() and self.getTransitDuration() and (planetaryRadius := self.getPlanetaryRadius()):
            return formulas.transitImpactParameter(self.radius, self.mass, planetaryRadius, self.period, self.transitDuration)
        return None

    def getOrbitalInclination(self):
        """Returns the orbital inclination axis of the planet with the lowest orbital period. 

        Returns:
            orbital inclination (float) - The angle between the star's and the planet's orbiting plane -- Units: `degrees`.
        """
        if self.radius and (semiMajorAxis := self.getSemiMajorAxis()) and (impactParameter := self.getImpactParameter()) is not None:
            return formulas.planetOrbitalInclination(self.radius, semiMajorAxis, impactParameter)
        return None
    
    def __calibrate(self, transitThreshold=1.5, timeStart=None, timeEnd=None, acceptanceLevel=ACCEPTANCE_LEVEL):
        """Identifies the correct phase and period.
        """
        if self.CALIBRATION_FLAG or transitThreshold < MINIMUM_TRANSIT_THRESHOLD: #Return if the transit threshold is too low or calibration has already occured.
            return
        #Finding the first TRANSIT_SAMPLES transits to test for the period.
        transitStart, transitPeak, transitEnd = self.transits.findTransitBounds(timeStart, transitThreshold=transitThreshold)
        transitSampleArray = [transitPeak]
        i = 0
        search_offset = CALIBRATION_SEARCH_OFFSET
        while transitPeak is not None and (i < TRANSIT_SAMPLES - 1 or (timeEnd is not None and transitPeak < timeEnd)): #Time is the limiting factor instead of the samples if it has been assigned a value.
            transitStart, transitPeak, transitEnd = self.transits.findTransitBounds(transitEnd + search_offset, timeEnd, transitThreshold=transitThreshold)
            if transitPeak is not None:
                transitSampleArray.append(transitPeak)
            i += 1
        #For each possible period (permutations of 2 transit samples), determine if the period is valid.
        for permutation in range(1,len(transitSampleArray)-1):
            for period, phase in ((period, transitSampleArray[i]) for i in range(len(transitSampleArray) - permutation - 1)
                                  #Permutation only checked if period is lower than the current period.
                                if (period := transitSampleArray[i + permutation] - transitSampleArray[i]) and MINIMUM_PERIOD < period and (self.period is None or self.period - search_offset)):
                #Initial probability of a valid period set to 50%.
                passed = 1
                total = 2
                #Terminate iterating if the probability exceeds or falls below the acceptance or rejection level.
                while REJECTION_LEVEL*total < passed < acceptanceLevel*total:
                    predictedTransit = phase + total*period
                    if self.transits.findTransit(predictedTransit - search_offset, predictedTransit + search_offset, transitThreshold=transitThreshold) is not None:
                        #If the transit threshold falls, more transits must pass.
                        passed += 1
                    total += 1
                #Conditions to check if new period is part of the old period by being a multiple of old period.
                if passed >= acceptanceLevel*total:
                    self.transitThreshold = transitThreshold
                    self.period = period
                    if self.phase is None or TRANSIT_CHARACTERISTICS_CLOSENESS_LEVEL < (self.period % period) < 1 - TRANSIT_CHARACTERISTICS_CLOSENESS_LEVEL or TRANSIT_CHARACTERISTICS_CLOSENESS_LEVEL < (self.phase - phase) % phase < 1 - TRANSIT_CHARACTERISTICS_CLOSENESS_LEVEL:
                        self.phase = phase
                    else:
                        self.period = period
            if self.period is not None: #No need to check for periods from larger permutations, as they would give larger periods.
                break
        if self.period is None and self.phase is None: #Performs calibration at a lower transit threshold if no valid period was identified.
            self.__calibrate(transitThreshold*TRANSIT_THRESHOLD_ITERATION_SCALING, acceptanceLevel=acceptanceLevel)
            if self.period is None and self.phase is None and acceptanceLevel > REJECTION_LEVEL:
                self.__calibrate(acceptanceLevel=acceptanceLevel - 0.2) #If no valid orbital periods found, sets the acceptance level to only 3 transits, then any transit pair.
        else: #Identifies if there is a lower valid period.
            self.__calibrate(transitThreshold*TRANSIT_THRESHOLD_ITERATION_SCALING, self.phase, self.phase + 2*self.period, acceptanceLevel=acceptanceLevel)
            self.CALIBRATION_FLAG = True

    def __sampleTransits(self):
        """Samples each transit by using a least squares sum formula to improve the orbital period and phase calculation.
        """
        if self.SAMPLING_FLAG:
            return
        self.__calibrate()
        nextTransitTimePredicted = self.phase + self.period
        lastTransitTime = self.transits.end
        nTransits = 1
        peakSum, weightedPeakSum = self.phase, 0
        nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum = 0, 0, 0
        searchOffset = self.period*self.transitThreshold*PERIOD_CALC_SEARCH_OFFSET
        recalibration = 2
        while nextTransitTimePredicted < lastTransitTime:
            nextTransitTimeFound = self.transits.findTransitPeak(nextTransitTimePredicted - searchOffset, nextTransitTimePredicted + searchOffset, transitThreshold=self.transitThreshold, searchMode=False)
            if nextTransitTimeFound is None:
                nSkippedTransits += 1
                skippedTransitsSum += nTransits
                skippedTransitsSquareSum += nTransits**2
                nextTransitTimePredicted += self.period
            else:
                peakSum += nextTransitTimeFound
                weightedPeakSum += nTransits*nextTransitTimeFound
                nextTransitTimePredicted = nextTransitTimeFound + self.period
            nTransits += 1
            if nTransits > recalibration:
                self.period = (nextTransitTimePredicted - self.phase)/nTransits
                recalibration *= RECALIBRATION_FACTOR

        self.period, self.phase = estimatePeriodicSignal(peakSum, weightedPeakSum, nTransits, nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum)
        self.SAMPLING_FLAG = True

    def __iter__(self):
        for dataHandler in LocalDataHandler():
            yield DataAnalyser(dataHandler=dataHandler)

class PhaseFoldedTransitModel():
    def __init__(self, phaseFoldedTimes, phaseFoldedFlux):
        """Creates a model for the phase-folded time-sorted transit data using polynomial interpolation.

        ----------
        Arguments:
            phaseFoldedTimes (arraylike) -- The phase-folded sorted time of the transit data.

            phaseFoldedFlux (arraylike) -- The flux sorted according to the phase-folded sorted time of the transit data.
        """
        #The phase-folded time-sorted transit data.
        self.phaseFoldedTimes, self.phaseFoldedFlux = phaseFoldedTimes, phaseFoldedFlux
        self.transitDetector = TransitDetector(phaseFoldedTimes, phaseFoldedFlux, searchMode=False)
        #Finds the transit bounds to create the interpolated model from.
        #Creates a polynomial to fit the phase folded transit.
        threshold = 1
        self.min, self.max = None, None
        while self.min is None or self.max is None:
            self.min, peak, self.max = self.transitDetector.findTransitBounds(0, searchMode=False, transitThreshold=threshold)
            threshold *= 0.75
        timeCFluxInterval = self.transitDetector[self.min:self.max]
        self.model = Polynomial.fit(*timeCFluxInterval, 4)
        #Finds the domain of the polynomial model.
        for root in sorted([x.real for x in self.model.roots() if np.isreal(x)]):
            if root > 0:
                self.max = root
                break
            else:
                self.min = root
        #Gets the coefficients of the polynomial model (used in evaluating the flux at a specified time in the __get_item__ function).
        self.coeffs = self.model.convert().coef
        self.peakFlux, self.peakTime = None, None

    def getPeak(self):
        """Returns the magnitude of the peak flux.

        Returns:
            peak flux (float) - The greatest change in flux when transiting -- Units: `Normalised stellar flux`.
        """
        if self.peakFlux is None:
            self.__initialisePeakValues()
        return -self.peakFlux
    
    def getPeakTime(self):
        """Returns the time at which the peak flux occured.

        Returns:
            time (float) -- Units: `Days`.
        """
        if self.peakTime is None:
            self.__initialisePeakValues()
        return self.peakTime

    def __initialisePeakValues(self):
        self.peakTime = min([x.real for x in polyroots(polyder(self.coeffs, 1)) if np.isreal(x) and self.min <= x <= self.max]) or self.min
        self.peakFlux = self[self.peakTime]
            
    def getData(self):
        """Returns the phase-folded time and the model's corresponding estimated flux values as a tuple of time and flux. 

        ----------
        Returns:
            tuple (phaseFoldedTimes, fluxModelEstimations):
              phaseFoldedTimes (arraylike) -- The phase-folded sorted time of the transit data.

              fluxModelEstimation (np.array) -- The flux evaluated from the model at the indexes corresponding to the phase-folded time.
        """
        return self.phaseFoldedTimes, np.fromiter(self, float)

    def __iter__(self):
        """Iterator of the flux evaluated from the model at the indexes corresponding to the phase-folded time.
        """
        for i in self.phaseFoldedTimes:
            yield self[i]

    def __getitem__(self, time):
        """Returns the model's estimated flux value of the phase folded transit data. 

        ----------
        Arguments:
            time (float) -- The time at which to evaluate the flux of the model.

        ----------
        Returns:
            flux (float) -- The flux at the specified time evaluated from the model.
        """
        if isinstance(time, slice):
            return np.fromiter((self[x] for x in np.arange(time.start, time.stop, time.step)), float)
        else:
            return polyval(time, self.coeffs) if self.min < time < self.max else .0