from math import floor
from statistics import variance
from data_handler import AbstractDataHandler, LocalDataHandler
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import Polynomial, polyval

SAMPLES = 200

def plot(*plots):
    plt.figure()
    plt.xlabel("Time")
    plt.ylabel("Flux")
    for data in plots:
        plt.plot(*data, '.', markersize=2)
    plt.show()

class TransitDetector():
    def __init__(self, times, flux, averageTime=False):
        self.times, self.flux = times, flux
        self.transitBound = None

        self.size = len(self.times)
        self.transitBound = None
        self.dt = min(self.times[1]-self.times[0], self.times[2]-self.times[1]) if averageTime else (self.times[-1] - self.times[0])/self.size
        self.standardStep = 1
        
        self.end = self.times[-1]
        self.start = self.times[0]

    def setStandardStep(self, time):
        self.standardStep = abs(floor(time/self.dt)) or 1

    def __findTime(self, time):
        return (self.times.searchsorted(time) or 1) - 1
    
    def __normaliseStep(self, step):
        return self.standardStep if step is None else floor(step*self.standardStep) or (self.standardStep if step >= 0 else -self.standardStep)

    def findTransitPeak(self, time=0, step=None):
        return self.times[self.__findTransitPeak(self.__findTime(time), self.__normaliseStep(step))]

    def __findTransitPeak(self, start=0, step=1):
        return self.__findTransitBounds(start, step)[1]
    
    def findTransitBounds(self, time=0, step=None):
        return tuple((self.times[x] for x in self.__findTransitBounds(self.__findTime(time), self.__normaliseStep(step))))

    def __findTransitBounds(self, start=0, step=1):
        """
        Returns the time index of a transit's start, peak and end.

        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            step (int) - The distance between samples for detecting a transit. When searching for a transit of a given length, the step 
            should be half the transit length.

        ----------
        Returns:
            tuple (start, peak, end):
                start (int) - The time index at which the transit is detected (flux is inside the transit bound consistently).
                peak (int) - The time index at which the transit reaches its maximum flux.
                end (int) - The time index at which the transit ends (flux is outside the transit bound consistently). 
        """
        start = self.__findTransitStart(start, step)
        fluxLow = self.flux[start]
        iLow = start
        end = self.size - 2 if step >= 0 else 0
        for i in range(start + 1, end, 1 if step >= 0 else -1):
            if (flux := self.flux[i]) < fluxLow:
                fluxLow, iLow = flux, i
            elif self.flux[i] > 0 and self.flux[i + 1] > 0 and self.flux[i + 2] > 0:
                return start, iLow, i
        raise Exception("Transit not found")

    def __findTransitStart(self, start=0, step=1):
        """
        Searches for the start of a transit. For a transit to be detected, the flux value must be inside the transit bound consistently.
        
        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            step (int) - The time index distance between samples for detecting a transit.

        ----------
        Returns:
            tuple (start, peak, end):
                start (int) - The time index at which the transit starts (flux is outside the transit bound consistently before this value).

                peak (int) - The time index at which the transit reaches its maximum flux.

                end (int) - The time index at which the transit ends (flux is outside the transit bound consistently after this value).
        """
        i = 0
        lb = self.getApproxTransitBound()
        end = self.size - 2 if step >= 0 else -1
        for i in range(start, end, step):
            if self.flux[i] < lb and self.flux[i + 1] < lb and self.flux[i + 2] < lb:
                break
        end = 1 if step >= 0 else self.size
        for i in range(i, end, -1 if step >= 0 else 1):
            if self.flux[i] > lb and self.flux[i - 1] > lb and self.flux[i - 2] > lb:
                return i
        raise Exception("Transit not found.")
            
    def getApproxTransitBound(self):
        """
        Returns an approximate flux value below which transits should occur:

        ----------
        Returns:
            flux value (float) -- Approx Transit Bound = Q1 - 2*IR, where Q1 is the lower quartile, IR the inequartile range.
        """
        if self.transitBound is None:
            sortedSamples = sorted(self.flux[:SAMPLES])
            self.transitBound = 3*sortedSamples[floor(0.25*SAMPLES)] - 2*sortedSamples[floor(0.75*SAMPLES)]
        return self.transitBound
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            start = self.__findTime(key.start)
            stop = self.__findTime(key.stop)
            return self.times[start:stop:key.step], self.flux[start:stop:key.step]
        else:
            i = self.__findTime(key)
            return self.times[i], self.flux[i]

class DataAnalyser():
    def __init__(self, dataID, dataHandler:AbstractDataHandler=LocalDataHandler):
        dataHandler = dataHandler(dataID)
        #Flux against Time data
        self.times, self.flux = dataHandler.getData()
        self.phaseFoldedTimes, self.phaseFoldedFlux = None, None
        self.transits = TransitDetector(self.times, self.flux)
        self.model = None
        #Stellar radius and mass
        self.radius, self.mass = dataHandler.getRadius(), dataHandler.getMass()


        self.size = len(self.times)
        self.period = None
        self.transitBound = None
        self.transitLength = None
        self.phase = None
        self.dt = min(self.times[1]-self.times[0], self.times[2]-self.times[1])

    def plot(self, plotType=""):
        match plotType:
            case "normal" | "standard" | "n" | "s":
                plot((self.times, self.flux))
            case "phase" | "phase folded" | "p":
                plot(self.getPhaseFoldedData())
            case "model" | "m":
                plot(self.getModel().getData())
            case "phase model" | "pm" | "p+":
                plot(self.getPhaseFoldedData(), self.getModel().getData())
        plt.show()

    def getData(self):
        return self.times, self.flux

    def getPhaseFoldedData(self):
        if self.phaseFoldedTimes is None:
            self.getOrbitalPeriod()
            self.phaseFoldedTimes = ((self.times - self.phase + self.period/2) % self.period) - self.period/2
            sort = np.argsort(self.phaseFoldedTimes)
            self.phaseFoldedTimes, self.phaseFoldedFlux = self.phaseFoldedTimes[sort], self.flux[sort]
        return self.phaseFoldedTimes, self.phaseFoldedFlux
    
    def getModel(self):
        if self.model is None:
            self.model = PhaseFoldedTransitModel(*self.getPhaseFoldedData())
        return self.model

    def getOrbitalPeriod(self):
        return self.period if self.period is not None else self.__calculateOrbitalPeriod()

    def getTransitLength(self):
        if self.transitLength is None:
            self.getPhase()
        return self.transitLength

    def getPhase(self):
        if self.phase is None:
            firstTransitStart, self.phase, firstTransitEnd = self.transits.findTransitBounds()
            self.transitLength = firstTransitEnd - firstTransitStart
        return self.phase
    
    def __calculateOrbitalPeriod(self):
        self.getPhase()
        self.period = self.phase - self.times[0]
        nextTransitTimePredicted = self.phase + self.period
        self.transits.setStandardStep(self.transitLength/2)
        lastTransit = self.transits.findTransitPeak(self.transits.end,-1)
        nTransitsStep, nTransits, skippedTransits = 1, 1, 0
        while nextTransitTimePredicted < lastTransit:
            nextTransitTimeFound = self.transits.findTransitPeak(nextTransitTimePredicted)
            if nTransits > 1:
                skippedTransits = round((nextTransitTimeFound-nextTransitTimePredicted)/self.period)
            self.period = (nextTransitTimeFound - self.phase)/(nTransits + skippedTransits)
            nTransitsStep *= 2
            nTransits += nTransitsStep + skippedTransits
            nextTransitTimePredicted = nextTransitTimeFound + self.period*(nTransitsStep-0.05)
        return (lastTransit-self.phase)/round((lastTransit-self.phase)/self.period)

class PhaseFoldedTransitModel():
    def __init__(self, phaseFoldedTimes, phaseFoldedFlux):
        """Creates a model for the phase-folded time-sorted transit data using polynomial interpolation.

        Arguments:
            phaseFoldedTimes (arraylike) -- The phase-folded sorted time of the transit data.

            phaseFoldedFlux (arraylike) -- The flux sorted according to the phase-folded sorted time of the transit data.
        """
        #The phase-folded time-sorted transit data.
        self.phaseFoldedTimes, self.phaseFoldedFlux = phaseFoldedTimes, phaseFoldedFlux
        self.transitDetector = TransitDetector(phaseFoldedTimes, phaseFoldedFlux, True)
        #Finds the transit bounds to create the interpolated model from.
        self.min, peak, self.max = self.transitDetector.findTransitBounds(0)
        #Creates a polynomial to fit the phase folded transit.
        self.model = Polynomial.fit(*self.transitDetector[self.min:self.max], 6)
        #Finds the domain of the polynomial model.
        for root in sorted([float(x) for x in self.model.roots() if np.isreal(x)]):
            if root > 0:
                self.max = root
                break
            else:
                self.min = root
        #Gets the coefficients of the polynomial model (used in evaluating the flux at a specified time in the __get_item__ function).
        self.coeffs = self.model.convert().coef

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
        return polyval(time, self.coeffs) if self.min < time < self.max else .0
    
    def chiSquaredFunc(self):
        N = len(self.phaseFoldedTimes)
        squaredValsGenerator = ((self.phaseFoldedFlux[i] - self[self.phaseFoldedTimes[i]])**2 for i in range(N))
        return (sum(squaredValsGenerator)/(N*variance(self.phaseFoldedFlux)))**0.5
