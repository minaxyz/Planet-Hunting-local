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
    t.out("Parameters")
    print(f"{period = }, {phase = }, {transitLength = }")
    t.totalOut()
    if plotType is not None:
        analyser.plot(plotType)

def iterTest():
    t = Timing(True, True)
    periods = {}
    print(f"{'Data ID':<15} | {'Period':<20} | {'Time'}")
    for d in DataAnalyser():
        periods[d.dataID] = d.getOrbitalPeriod()
        print(f"{d.dataID:<15} | {periods[d.dataID]:<20} | ", end=" ")
        t.out()
    t.totalOut()
    
#KIC002571238 period = 9.286958783276173
#kplr002715135 period = 5.74771
#timedTest("kplr002715135", "pm")
iterTest()