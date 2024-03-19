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
    t.out("Parameters")
    print(f"{period = }, {phase = }, {transitLength = }, {peak = }, {threshold = }")
    print(analyser.transits.getAnomalousRegions())
    t.totalOut()
    if plotType is not None:
        analyser.plot(plotType)

def iterTest():
    t = Timing(True, True)
    failed = []
    print(f"{'Data ID':<15} | {'Period':<20} | {'Planetary Radius':<20} | {'Time'}")
    for d in DataAnalyser():
        try:
            print(f"{d.dataID:<15} | {d.getOrbitalPeriod():<20} | {d.getPlanetaryRadius():<20} | ", end=" ")
            t.out()
        except Exception:
            failed.append(d.dataID)
    if failed:
        print("Failed:")
        for dataID in failed:
            print(dataID)
    else:
        print("No fails.")
    t.totalOut()
    
#KIC002571238 period = 9.286958783276173
#kplr002715135 period = 5.74771
#KPLR009266431 period = 18.3963
#timedTest("KPLR009266431", "n")
iterTest()