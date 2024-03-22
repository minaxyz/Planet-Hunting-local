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
TRANSIT_THRESHOLD_ITERATION_SCALING = 0.75 #The scaling applied to the threshold after each iteration. Must be between 0 and 1.
MINIMUM_TRANSIT_THRESHOLD = 0.5 #Determines the maximum number of recursion of calibration (max recursion depth = ceil(log_TRANSIT_THRESHOLD_ITERATION_SCALING(MINIMUM_TRANSIT_THRESHOLD))).
"""Determines the strictness in validating a transit (Higher -> Higher strictness), must be between 0 and 1, higher than ACCEPTANCE_LEVEL.
Currently, accepts a transit after 3 passes and 0 fails, or accepted if above 80% passes. (#Passes no fails to accept = (1/1 - ACCEPTANCE_LEVEL) - 2))"""
ACCEPTANCE_LEVEL = 0.8
"""Determines the strictness in rejecting a transit (Lower -> Higher strictness), must be between 0 and 1, lower than REJECTION_LEVEL.
Currently, rejects a transit after 0 passes and 3 fails, or rejected if below 20% passes. (#Fails no passes to reject = (1/REJECTION_LEVEL) - 2"""
REJECTION_LEVEL = 0.2
##Orbital Period Calculation Constants
SEARCH_OFFSET = 0.1 #Determines the interval to search for the transit.

def plot(*plots):
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("Flux")
    for data in plots:
        plt.plot(*data, '.', markersize=2)
    plt.show()

def histogram(flux):
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("Proportion")
    plt.hist(flux, 100, weights=np.full(len(flux),1/len(flux),dtype=float))
    plt.show()

def phaseFold(times, flux, period, phase):
    phaseFoldedTimes = ((times - phase + period/2) % period) - period/2
    sort = np.argsort(phaseFoldedTimes)
    return phaseFoldedTimes[sort], flux[sort]

class TransitDetector():
    def __init__(self, times, flux, searchMode=True):
        self.times = times
        self.size = len(self.times)
        self.dt = (self.times[-1] - self.times[0])/self.size
        #Applies a convolution to the flux which reduces noise and accentuates transits.
        self.uniformConvolutionFactor = floor((0.5 if searchMode else 0.05)*MINIMUM_PERIOD/self.dt)
        self.convolvedFlux = np.convolve(flux, np.full(self.uniformConvolutionFactor, 1/self.uniformConvolutionFactor, dtype=float),'same')

        self.anomalousRegions = [] #Even indexes represent the start of an anomalous region, and its next odd index represents the end of the anomaly.
        self.medianFlux = None #The median flux, calculated from the first FLUX_SAMPLES flux values.
        self.transitThreshold = None #The threshold below which transits should occur. Implementation details in self.getTransitThreshold.
        
        self.end = self.times[-1]
        self.start = self.times[0]

    def __findTime(self, time):
        return min(self.times.searchsorted(time), self.size - 1)

    def findTransit(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None):
        self.getTransitThreshold()
        i = self.__findTransit(
            self.__findTime(timeStart) if timeStart is not None else (self.size if reverse else 0),
            self.__findTime(timeEnd) if timeEnd is not None else (-1 if reverse else self.size), 
            reverse, 
            self.transitThreshold if transitThreshold is None else self.medianFlux + transitThreshold*(self.transitThreshold - self.medianFlux))
        return i and self.times[i]

    def findTransitPeak(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None, searchMode=True):
        return self.findTransitBounds(timeStart, timeEnd, reverse, transitThreshold, searchMode)[1]
    
    def findTransitBounds(self, timeStart=None, timeEnd=None, reverse=False, transitThreshold=None, searchMode=True):
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
                return leftTransitBounds if 2*midPoint > leftTransitBounds[1] + rightTransitBounds[1] else rightTransitBounds
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

                None is returned if a transit is not found.
        """
        START = start
        END = end
        #Ensuring that the transit is not anomalous.
        while (start := self.__findTransit(start, end, reverse, transitThreshold)) != (nregion := self.__findNormalRegion(start, reverse)):
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
        self.addAnomalousRegion(*((anomalyEnd, firstAnomaly) if reverse else (firstAnomaly, anomalyEnd)))
        return anomalyEnd
    
    def addAnomalousRegion(self, start, end):
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
        return [[self.times[x] for x in self.anomalousRegions[i:i+2]] for i in range(0,len(self.anomalousRegions),2)]

    def getData(self):
        return self.times, self.convolvedFlux

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = self.__findTime(key.start)
            stop = self.__findTime(key.stop)
            return self.times[start:stop], self.convolvedFlux[start:stop]
        else:
            i = self.__findTime(key)
            return self.times[i], self.convolvedFlux[i]

class DataAnalyser():

    def __init__(self, dataID=None, dataHandler:AbstractDataHandler=LocalDataHandler):
        self.dataHandler = dataHandler(dataID) if dataID else dataHandler if not inspect.isclass(dataHandler) else None
        self.dataID = dataID or self.dataHandler and self.dataHandler.dataID
        #Flux against Time data
        self.times, self.flux = self.dataHandler.getData() if self.dataHandler is not None else (None, None)


        self.phaseFoldedTimes, self.phaseFoldedFlux = None, None
        self.transits = TransitDetector(self.times, self.flux) if self.dataHandler is not None else None
        self.model = None
        #Stellar radius and mass
        self.radius, self.mass = (self.dataHandler.getRadius(), self.dataHandler.getMass()) if self.dataHandler is not None else (None, None)

        self.CALIBRATION_FLAG = False
        self.ACCURATE_PERIOD_FLAG = False
        self.size = len(self.times) if self.dataHandler is not None else None
        self.period = None
        self.transitThreshold = None
        self.transitLength = None
        self.phase = None

    def plot(self, plotType=""):
        match plotType:
            case "normal" | "standard" | "n" | "s" | "":
                plot(self.getData())
            case "phase" | "phase folded" | "p":
                plot(self.getPhaseFoldedData())
            case "convolved" | "convolution" | "con" | "c":
                plot((self.times, self.transits.convolvedFlux))
            case "model" | "m":
                plot(self.getModel().getData())
            case "phase model" | "pm" | "p+":
                plot(self.getPhaseFoldedData(), self.getModel().getData())
            case "convolved phase" | "cp":
                plot(self.getModel().transitDetector.getData())
            case "convolved phase model" | "cpm":
                plot(self.getModel().transitDetector.getData(), self.getModel().getData())
            case "histogram" | "hist" | "h":
                histogram(self.flux)
            case "convolved histogram" | "con hist" | "ch":
                histogram(self.transits.convolvedFlux)
            case _:
                raise Exception("Invalid Plot Type: Plot not recognised.\nPlot Type Options: 'standard', 'phase folded', 'model', 'phase model'.")
        plt.show()

    def getData(self):
        return self.times, self.flux

    def getPhaseFoldedData(self):
        if self.phaseFoldedTimes is None:
            self.getOrbitalPeriod()
            self.phaseFoldedTimes, self.phaseFoldedFlux = phaseFold(self.times, self.flux, self.period, self.phase)
        return self.phaseFoldedTimes, self.phaseFoldedFlux
    
    def getModel(self):
        if self.model is None:
            self.model = PhaseFoldedTransitModel(*self.getPhaseFoldedData())
            self.transitLength = self.model.max - self.model.min
        return self.model
    
    
    def getOrbitalPeriod(self):
        if not self.ACCURATE_PERIOD_FLAG:
            self.__calculateOrbitalPeriod()
        return self.period

    def getTransitLength(self):
        if self.transitLength is None:
            start, peak, end = self.transits.findTransitBounds(self.getPhase(), transitThreshold=self.transitThreshold)
            if start and end:
                self.transitLength = max(end - start - 2*self.transits.dt*self.transits.uniformConvolutionFactor, self.transits.dt)
        return self.transitLength

    def getPhase(self):
        if self.phase is None:
            self.__calibrate()
        return self.phase
    
    def getPlanetaryRadius(self):
        return formulas.planetaryRadius(self.radius, self.getModel().getPeak())
    
    def getSemiMajorAxis(self):
        return formulas.semiMajorAxis(self.mass, self.getOrbitalPeriod())

    def getImpactParameter(self):
        return formulas.transitImpactParameter(self.radius, self.getPlanetaryRadius(), self.getOrbitalPeriod(), self.getTransitLength())

    def getOrbitalInclination(self):
        return formulas.planetOrbitalInclination(self.radius, self.getPlanetaryRadius(), self.mass, self.getOrbitalPeriod(), self.getTransitLength())
    
    def __calibrate(self, transitThreshold=1, timeStart=None, timeEnd=None):
        """Initialises the period.
        """
        if transitThreshold < MINIMUM_TRANSIT_THRESHOLD or self.CALIBRATION_FLAG: #Transit threshold is too low or calibration has already occured.
            return
        #Finding the first TRANSIT_SAMPLES transits to test for the period.
        transitStart, transitPeak, transitEnd = self.transits.findTransitBounds(timeStart, transitThreshold=transitThreshold)
        transitSampleArray = [transitPeak]
        i = 0
        search_offset = 2*self.transits.dt
        while i < TRANSIT_SAMPLES - 1 and timeEnd is not None: #Time is the limiting factor instead of the samples if it has been assigned a value.
            transitStart, transitPeak, transitEnd = self.transits.findTransitBounds(transitEnd + search_offset, timeEnd, transitThreshold=transitThreshold)
            if transitPeak is not None:
                transitSampleArray.append(transitPeak)
            else:
                break
            i += 1
        #For each possible period (permutations of 2 transit samples), determine if the period is valid.
        for permutation in range(1,len(transitSampleArray)-1):
            for period, phase in ((period, transitSampleArray[i]) for i in range(len(transitSampleArray) - permutation - 1)
                                  #Permutation only checked if period is lower than the current period.
                                if (period := transitSampleArray[i + permutation] - transitSampleArray[i]) and self.period is None or period > MINIMUM_PERIOD and period < self.period - search_offset):
                #Initial probability of a valid period set to 50%.
                passed = 1
                total = 2
                #Terminate iterating if the probability exceeds or falls below the acceptance or rejection level.
                while REJECTION_LEVEL*total < passed < ACCEPTANCE_LEVEL*total:
                    predictedTransit = phase + total*period
                    if (transit := self.transits.findTransit(predictedTransit - search_offset, predictedTransit + search_offset, transitThreshold=transitThreshold)) is not None:
                        #If the transit threshold falls, more transits must pass.
                        passed += 1
                    total += 1
                if passed >= ACCEPTANCE_LEVEL*total:
                    self.transitThreshold = transitThreshold
                    self.period = period
                    #Conditions to determing whether the phase of the old period if part of the new period. 
                    if self.phase is None and abs(self.period-period)/self.period < 0.0125:
                        self.phase = phase
                    else:
                        n, diff = divmod(abs(self.phase-phase),period)
                        if diff > n*0.05:
                            self.phase = phase
            if self.period is not None: #No need to check for periods from larger permutations, as they would give larger periods.
                break
        if self.period is None and self.phase is None: #Performs calibration at a lower transit threshold if no valid period was identified.
            self.__calibrate(transitThreshold*TRANSIT_THRESHOLD_ITERATION_SCALING)
        else: #Identifies if there is a lower valid period.
            self.__calibrate(transitThreshold*TRANSIT_THRESHOLD_ITERATION_SCALING, self.phase, self.phase + 2*self.period)
            self.CALIBRATION_FLAG = True


    def __calculateOrbitalPeriod(self):
        """Uses a least squares sum method to calculate the orbital period, and improves the estimation of the phase.
        """
        self.__calibrate()
        nextTransitTimePredicted = self.phase + self.period
        lastTransit = self.transits.findTransitPeak(self.transits.end, reverse=True, transitThreshold=self.transitThreshold)
        nTransits = 1
        peakSum, weightedPeakSum = self.phase, 0
        nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum = 0, 0, 0
        backtrack = max(self.period*self.transitThreshold*SEARCH_OFFSET, self.getTransitLength())
        recalibration = 2
        while nextTransitTimePredicted < lastTransit:
            nextTransitTimeFound = self.transits.findTransitPeak(nextTransitTimePredicted - backtrack, nextTransitTimePredicted + backtrack, transitThreshold=self.transitThreshold, searchMode=False)
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
            if nTransits == recalibration:
                self.period = (nextTransitTimePredicted - self.phase)/nTransits
                recalibration *= 2

        self.period, self.phase = estimatePeriodicSignal(peakSum, weightedPeakSum, nTransits, nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum)
        self.ACCURATE_PERIOD_FLAG = True

    def __iter__(self):
        for dataHandler in LocalDataHandler():
            yield DataAnalyser(dataHandler=dataHandler)

class PhaseFoldedTransitModel():
    def __init__(self, phaseFoldedTimes, phaseFoldedFlux):
        """Creates a model for the phase-folded time-sorted transit data using polynomial interpolation.

        Arguments:
            phaseFoldedTimes (arraylike) -- The phase-folded sorted time of the transit data.

            phaseFoldedFlux (arraylike) -- The flux sorted according to the phase-folded sorted time of the transit data.
        """
        #The phase-folded time-sorted transit data.
        self.phaseFoldedTimes, self.phaseFoldedFlux = phaseFoldedTimes, phaseFoldedFlux
        self.transitDetector = TransitDetector(phaseFoldedTimes, phaseFoldedFlux, searchMode=False)
        #Finds the transit bounds to create the interpolated model from.
        self.min, peak, self.max = self.transitDetector.findTransitBounds(0, searchMode=False)
        #Creates a polynomial to fit the phase folded transit.
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
        if self.peakFlux is None:
            self.__initialisePeakValues()
        return self.peakFlux
    
    def getPeakTime(self):
        if self.peakTime is None:
            self.__initialisePeakValues()
        return self.peakTime

    def __initialisePeakValues(self):
        peakTimesArray = [x.real for x in polyroots(polyder(self.coeffs, 1)) if np.isreal(x)]
        self.peakTime = peakTimesArray[0]
        for peakTime in peakTimesArray[1:]:
            if abs(peakTime) < abs(self.peakTime):
                self.peakTime = peakTime
        self.peakFlux = self[self.peakTime]
            
    def getData(self):
        """Returns the phase-folded time and the model's corresponding estimated flux values as a tuple of time and flux. 

        Returns:
            tuple (phaseFoldedTimes, fluxModelEstimations)
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

        Arguments:
            time (float) -- The time at which to evaluate the flux of the model.

        Returns:
            flux (float) -- The flux at the specified time evaluated from the model.
        """
        if isinstance(time, slice):
            return np.fromiter((self[x] for x in np.arange(time.start, time.stop, time.step)), float)
        else:
            return polyval(time, self.coeffs) if self.min < time < self.max else .0