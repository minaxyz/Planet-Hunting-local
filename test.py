from data_handler import LocalDataHandler
from data_analyser import DataAnalyser
from time import perf_counter
import matplotlib.pyplot as plt
import formulas

class Timing():
    def __init__(self, sinceLastOut=True, ms=True):
        self.s = perf_counter()
        self.last = self.s
        self.sinceLastOut = sinceLastOut
        self.unit = 1000 if ms else 1
        self.unitPrefix = 'm'*ms + 's'

    def out(self,label=None):
        print(f"{label + ':' if label else ''} {((e := perf_counter()) - self.last)*self.unit} {self.unitPrefix}")
        if self.sinceLastOut:
            self.last = e

    def totalOut(self, label=""):
        if not self.sinceLastOut:
            self.out(label)
        print(f"Total: {(perf_counter() - self.s)*self.unit} {self.unitPrefix}")
    
    def skip(self):
        self.last = perf_counter()

def timedTest(dataID, plotType=None):
    print(f"{dataID} results:")
    t = Timing(True, True)
    analyser = DataAnalyser(dataID)
    t.out("Initialisation")
    transitLength = analyser.getTransitLength()
    phase = analyser.getPhase()
    period = analyser.getOrbitalPeriod()
    peak = analyser.getModel().getPeak()
    threshold = analyser.transits.getTransitThreshold()
    orbitalInclination = analyser.getOrbitalInclination()
    impactParameter = analyser.getImpactParameter()
    semiMajorAxis = analyser.getSemiMajorAxis()
    t.out("Parameters")
    print(f"{period = }, {phase = }, {transitLength = }, {peak = }, {threshold = }, {orbitalInclination = }, {semiMajorAxis = }, {impactParameter = }")
    print(f"Anomalous Flux Regions: {analyser.transits.getAnomalousRegions()}")
    t.totalOut()
    if plotType is not None:
        analyser.plot(plotType)

def iterTest():
    t = Timing(True, True)
    failed = []
    i = 0
    print(f"{'N':<6} | {'Data ID':<15} | {'Period':<20} | {'Planetary Radius':<20} | {'Semi Major Axis':<20} | {'Impact Parameter':<20} | {'Orbital Inclination':<20} | {'Time'}")
    for d in DataAnalyser():
        i += 1
        try:
            d.getModel()
            print(f"{i:<6} | {d.dataID:<15} | {d.getOrbitalPeriod():<20} | {d.getPlanetaryRadius():<20} | {d.getSemiMajorAxis():<20} | {d.getImpactParameter():<20} | {d.getOrbitalInclination():<20} |", end=" ")
            t.out()
        except Exception:
            failed.append(d.dataID)
            t.skip()
    if failed:
        print(f"Failed: {len(failed)}/{i}")
        for dataID in failed:
            print(dataID)
    else:
        print(f"No fails. 0/{i}")
    t.totalOut()

#DataAnalyser("kplr005617854").plot('c')
#TODO: Revisit kplr005617854
timedTest("KIC002571238", "p")
#iterTest()