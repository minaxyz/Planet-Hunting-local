import matplotlib.pyplot as plt
import formulas
import numpy as np
import inspect
from numpy.polynomial.polynomial import Polynomial, polyval, polyder, polyroots
from math import floor

from data_handler import AbstractDataHandler, LocalDataHandler
from utils import estimatePeriodicSignal

SAMPLES = 300

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
    plt.ylabel("Frequency")
    plt.hist(flux, 100)
    plt.show()

class TransitDetector():
    def __init__(self, times, flux, averageTime=False):
        self.times, self.flux = times, flux

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
        end = self.size - 2
        i = start
        for i in range(start + 1, end, 1):
            if (flux := self.flux[i]) < fluxLow:
                fluxLow, iLow = flux, i
            elif self.flux[i] > 0 and self.flux[i + 1] > 0 and self.flux[i + 2] > 0:
                return start, iLow, i
        return start, iLow, i

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
        if start > self.size - 3:
            start = self.size - 3
        end = self.size - 2 if step >= 0 else -1
        for i in range(start, end, step):
            if self.flux[i] < lb and self.flux[i + 1] < lb and self.flux[i + 2] < lb:
                break
        end = -1 if step >= 0 else self.size
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
            self.transitBound = -sortedSamples[floor(0.99*SAMPLES)]
        return self.transitBound
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            start = self.__findTime(key.start)
            stop = self.__findTime(key.stop)
            return self.times[start:stop], self.flux[start:stop]
        else:
            i = self.__findTime(key)
            return self.times[i], self.flux[i]

class DataAnalyser():
    def __init__(self, dataID=None, dataHandler:AbstractDataHandler=LocalDataHandler):
        if dataID is not None or not inspect.isclass(dataHandler):
            self.load(dataID, dataHandler)
        else:
            self.times, self.flux, self.radius, self.mass, self.period, self.phase, self.transitLenth, self.transits, self.model, self.phaseFoldedTimes, self.phaseFoldedFlux = None, None, None, None, None, None, None, None, None, None, None
        

    def load(self, dataID=None, dataHandler:AbstractDataHandler=LocalDataHandler):
        self.dataHandler = dataHandler(dataID) if inspect.isclass(dataHandler) else dataHandler
        #Flux against Time data
        self.times, self.flux = self.dataHandler.getData()

        self.phaseFoldedTimes, self.phaseFoldedFlux = None, None
        self.transits = TransitDetector(self.times, self.flux)
        self.model = None
        #Stellar radius and mass
        self.radius, self.mass = self.dataHandler.getRadius(), self.dataHandler.getMass()


        self.size = len(self.times)
        self.period = None
        self.transitLenth = None
        self.phase = None

    def plot(self, plotType=""):
        match plotType:
            case "normal" | "standard" | "n" | "s" | "":
                plot(self.getData())
            case "phase" | "phase folded" | "p":
                plot(self.getPhaseFoldedData())
            case "model" | "m":
                plot(self.getModel().getData())
            case "phase model" | "pm" | "p+":
                plot(self.getPhaseFoldedData(), self.getModel().getData())
            case "histogram" | "hist" | "h":
                histogram(self.flux)
            case _:
                raise Exception("Invalid Plot Type: Plot not recognised.\nPlot Type Options: 'standard', 'phase folded', 'model', 'phase model'.")
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
            self.transitLenth = self.model.max - self.model.min
        return self.model
    
    def getOrbitalPeriod(self):
        if self.period is None:
            self.__calculateOrbitalPeriod()
        return self.period

    def getTransitLength(self):
        if self.transitLength is None:
            self.getPhase()
        return self.transitLength

    def getPhase(self):
        if self.phase is None:
            firstTransitStart, self.phase, firstTransitEnd = self.transits.findTransitBounds()
            self.transitLenth = firstTransitEnd - firstTransitStart
        return self.phase
    
    def getPlanetaryRadius(self):
        return formulas.PlanetaryRadius(self.mass, self.getModel().getPeak())


    def __calculateOrbitalPeriod(self):
        """Uses a least squares sum method to calculate the orbital period, and improves the estimation of the phase.
        """
        self.getPhase()
        self.transits.setStandardStep(self.transitLenth/6)
        self.period = self.phase - self.times[0]
        nextTransitTimePredicted = self.phase + self.period
        lastTransit = self.transits.findTransitBounds(self.transits.end,-1)[0]
        backtrack = -0.05*self.period
        nTransits = 1
        peakSum, weightedPeakSum = self.phase, 0
        nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum = 0, 0, 0
        while nextTransitTimePredicted < lastTransit:

            nextTransitTimeFound = self.transits.findTransitPeak(nextTransitTimePredicted + backtrack)
            if (currentSkippedTransits := round((nextTransitTimeFound-nextTransitTimePredicted)/self.period)) and nTransits > 1:
                newNTransits = nTransits + currentSkippedTransits
                nSkippedTransits += currentSkippedTransits
                skippedTransitsSum += sum((x for x in range(nTransits,newNTransits)))
                skippedTransitsSquareSum += sum((x*x for x in range(nTransits,newNTransits)))
                nTransits = newNTransits
            peakSum += nextTransitTimeFound
            weightedPeakSum += nTransits*nextTransitTimeFound
            self.period = (nextTransitTimeFound - self.phase)/nTransits
            nextTransitTimePredicted = nextTransitTimeFound + self.period
            nTransits += 1
        
        self.period, self.phase = estimatePeriodicSignal(peakSum, weightedPeakSum, nTransits, nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum)

    def __iter__(self):
        for handler in [x for x in LocalDataHandler() if x.dataID.startswith('KIC')]:
            yield DataAnalyser(dataHandler=handler)

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
        var = self.transitDetector[self.min:self.max]
        self.model = Polynomial.fit(*var, 4)
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
