import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import Polynomial, polyval
from math import floor

from data_handler import AbstractDataHandler, LocalDataHandler
from utils import estimatePeriodicSignal

SAMPLES = 300
SIGNIFICANCE_LEVEL = 0.05
MINIMUM_PERIOD = 0.241842217

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

class TransitDetector():
    def __init__(self, times, flux, averageTime=False):
        self.times = times
        self.size = len(self.times)
        self.dt = min(self.times[1]-self.times[0], self.times[2]-self.times[1]) if averageTime else (self.times[-1] - self.times[0])/self.size
        #Applies a convolution to the flux which reduces noise and accentuates transits.
        uniformConvolutionFactor = floor(0.95*MINIMUM_PERIOD/self.dt) or 1
        self.convolutedFlux = np.convolve(flux, np.full(uniformConvolutionFactor, 1/uniformConvolutionFactor, dtype=float),'same')

        self.transitBound = None
        self.standardStep = 1
        
        self.end = self.times[-1]
        self.start = self.times[0]

    def __findTime(self, time):
        return min(self.times.searchsorted(time), self.size - 1)

    def findTransitPeak(self, time=0, reverse=False):
        return self.times[self.__findTransitPeak(self.__findTime(time), reverse)]

    def __findTransitPeak(self, start, reverse):
        return self.__findTransitBounds(start, reverse)[1]
    
    def findTransitBounds(self, time=0, reverse=False):
        return tuple((self.times[x] for x in self.__findTransitBounds(self.__findTime(time), reverse)))

    def __findTransitBounds(self, start, reverse):
        """
        Returns the time index of a transit's start, peak and end.

        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            revese (bool) - Reverses the direction of the search if True, by default the timestep is the same as the time array.

        ----------
        Returns:
            tuple (start, peak, end):
                start (int) - The time index at which the transit is detected.
                
                peak (int) - The time index at which the transit reaches its maximum flux.
                
                end (int) - The time index at which the transit ends. 
        """
        start = self.__findTransit(start, reverse)
        for i in range(start, -1, -1):
            if self.convolutedFlux[i] > self.transitBound:
                break
        fluxLow = self.convolutedFlux[start]
        iLow = start
        i = start
        for i in range(start + 1, self.size, -1 if reverse else 1):
            if (flux := self.convolutedFlux[i]) < fluxLow:
                fluxLow, iLow = flux, i
            elif self.convolutedFlux[i] > self.transitBound:
                return start, iLow, i
        return start, iLow, i

    def __findTransit(self, start, reverse):
        """
        Searches for the start of a transit. For a transit to be detected, the convolved flux value must be below the transit bound.
        
        ----------
        Parameters:
            start (int) - The time index the searching of a transit begins from.

            revese (bool) - Reverses the direction of the search if True, by default the timestep is the same as the time array.

        ----------
        Returns:
            transit bound (int) - The time index at which the transit is detected.
        """
        i = 0
        self.getApproxTransitBound()
        if start > self.size - 3:
            start = self.size - 3
        start = min(max(start, 0), self.size - 3)
        end = -1 if reverse else self.size
        for i in range(start, end, -1 if reverse else 1):
            if self.convolutedFlux[i] < self.transitBound:
                return i
        raise Exception("Transit not found.")
            
    def getApproxTransitBound(self):
        """
        Returns an approximate flux value below which transits should occur.

        ----------
        Returns:
            flux value (float) -- Defined as `Median convolved flux value - 2.5(SIGNIFICANCE_LEVEL from the highest convolved flux value)`
            from a sample of the first SAMPLE values.
        """
        if self.transitBound is None:
            sortedSamples = sorted(self.convolutedFlux[:SAMPLES])
            self.transitBound = sortedSamples[floor(0.5*SAMPLES)] - 2.5*sortedSamples[floor((1 - SIGNIFICANCE_LEVEL)*SAMPLES)]
        return self.transitBound
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            start = self.__findTime(key.start)
            stop = self.__findTime(key.stop)
            return self.times[start:stop], self.convolutedFlux[start:stop]
        else:
            i = self.__findTime(key)
            return self.times[i], self.convolutedFlux[i]

class DataAnalyser():
    def __init__(self, dataID, dataHandler:AbstractDataHandler=LocalDataHandler):
        self.dataHandler = dataHandler(dataID)
        #Flux against Time data
        self.times, self.flux = self.dataHandler.getData()
        self.phaseFoldedTimes, self.phaseFoldedFlux = None, None
        self.transits = TransitDetector(self.times, self.flux)
        self.model = None
        #Stellar radius and mass
        self.radius, self.mass = self.dataHandler.getRadius(), self.dataHandler.getMass()


        self.size = len(self.times)
        self.period = None
        self.transitBound = None
        self.transitLenth = None
        self.phase = None
        self.dt = min(self.times[1]-self.times[0], self.times[2]-self.times[1])

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
            case "convolved" | "convolution" | "con" | "c":
                plot((self.times, self.transits.convolutedFlux))
            case "histogram" | "hist" | "h":
                histogram(self.flux)
            case "convolved histogram" | "con hist" | "ch":
                histogram(self.transits.convolutedFlux)
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
    
    def __calculateOrbitalPeriod(self):
        """Uses a least squares sum method to calculate the orbital period, and improves the estimation of the phase.
        """
        self.getPhase()
        self.period = self.transits.findTransitPeak(self.phase + self.transitLenth) - self.phase
        nextTransitTimePredicted = self.phase + self.period
        lastTransit = self.transits.findTransitBounds(self.transits.end, True)[0]
        nTransits = 1
        peakSum, weightedPeakSum = self.phase, 0
        nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum = 0, 0, 0
        while nextTransitTimePredicted < lastTransit:
            nextTransitTimeFound = self.transits.findTransitPeak(nextTransitTimePredicted)
            if (currentSkippedTransits := round((nextTransitTimeFound-nextTransitTimePredicted)/self.period)) and nTransits > 1:
                newNTransits = nTransits + currentSkippedTransits
                nSkippedTransits += currentSkippedTransits
                skippedTransitsSum += sum((x for x in range(nTransits,newNTransits)))
                skippedTransitsSquareSum += sum((x*x for x in range(nTransits,newNTransits)))
                nTransits = newNTransits
            peakSum += nextTransitTimeFound
            weightedPeakSum += nTransits*nextTransitTimeFound
            self.period = (nextTransitTimeFound - self.phase)/nTransits
            nextTransitTimePredicted = nextTransitTimeFound + (1 - SIGNIFICANCE_LEVEL)*self.period
            nTransits += 1
        
        self.period, self.phase = estimatePeriodicSignal(peakSum, weightedPeakSum, nTransits, nSkippedTransits, skippedTransitsSum, skippedTransitsSquareSum)

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
