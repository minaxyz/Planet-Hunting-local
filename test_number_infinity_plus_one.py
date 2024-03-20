from data_analyser import TransitDetector, DataAnalyser
from data_handler import AbstractDataHandler, LocalDataHandler
import matplotlib.pyplot as plt

DataTitles = ['KIC002571238','KIC005881688','KIC006922244','KIC007950644','KIC008359498','KIC009631995','KIC010418224','KIC011853905']
tranistlengths = []

for i in range(len(DataTitles)):
    planet = DataAnalyser(DataTitles[i])
    tranistlengths.append(planet.getTransitLength())

print(tranistlengths)
plt.hist(tranistlengths)
plt.xlabel('Transit lengths')
plt.ylabel('Number of planets')
plt.show()