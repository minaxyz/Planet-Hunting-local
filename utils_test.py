from utils import *
from time import perf_counter
import random
import numpy as np

def benchmarkArrayFunc(size, n, funcs):
    testList = np.asarray([x + random.random() for x in range(size)])
    elements = [0 for x in range(n)]
    resultList = [[] for x in range(len(funcs))]
    tSum = [0 for x in funcs]
    for iteration in range(n):
        element = random.randint(0,size-1) + random.random()
        elements[iteration] = element
        for i, func in enumerate(funcs):
            t1 = perf_counter()
            result = func(testList, element)
            t2 = perf_counter()
            resultList[i].append(result)
            tSum[i] += t2 - t1
    for i, func in enumerate(funcs):
        print(f"{func.__name__}: {tSum[i]*1000/n}ms")

def testPeriodicSignalEstimator(period, phase, size, absSpread=0, missingPercent=0):
    array = [period*x + phase + 2*absSpread*random.random() - absSpread for x in range(size)]
    missingPoints = random.sample([x for x in range(size)], floor(size*missingPercent))
    xSum = sum([array[i] for i in range(size) if i not in missingPoints])
    xiSum = sum((i*array[i] for i in range(size) if i not in missingPoints))
    N = size
    s1 = len(missingPoints)
    s2 = sum(missingPoints)
    s3 = sum([x*x for x in missingPoints])
    periodE, phaseE = estimatePeriodicSignal(xSum, xiSum, N, s1, s2, s3)
    print(f"Period - Actual: {period}, Estimate: {periodE}\nPhase - Actual: {phase}, Estimate: {phaseE}")

#benchmarkArrayFunc(65000, 10000, [binarySearch, np.searchsorted, interpolationSearch])
testPeriodicSignalEstimator(30, 100, 100, 1, 0.2)