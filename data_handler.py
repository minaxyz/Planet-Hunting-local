import formulas
from abc import ABC, abstractmethod
from PyPDF2 import PdfReader
import numpy as np
import os

class AbstractDataHandler(ABC):
    def __init__(self, dataID):
        pass

    @abstractmethod
    def getData(self):
        return None

class LocalDataHandler(AbstractDataHandler):
    systemsDirectoryDict = None
    stellarMassRadiusDict = None

    def __init__(self, dataID=None):
        if self.systemsDirectoryDict is None:
            self.__initialiseSystemsDirectoryDict()
        if self.stellarMassRadiusDict is None:
            self.__initialiseStellarMassAndRadius()
        if dataID is None:
            self.dataID, self.times, self.flux, self.radius, self.mass = None, None, None, None, None
        else:
            self.load(dataID)

    def __initialiseSystemsDirectoryDict(self):
        self.systemsDirectoryDict = {file.split('.')[0].upper():"Data/" + file for file in os.listdir("Data") if file.endswith((".tbl",".dat")) and "phaseFold" not in file and "TIC" not in file
                                } | {file.split('_')[0].upper():"NewData/" + file for file in os.listdir("NewData") if file.endswith(".tbl")}
    
    def __initialiseStellarMassAndRadius(self):
        f = PdfReader(open("Data/Stellar_Mass_Radius.pdf", 'rb'))
        tableData = f.pages[0].extract_text().split()[9:]
        self.stellarMassRadiusDict = {tableData[i]:(float(tableData[i+1]),float(tableData[i+2]))  for i in range(0,len(tableData),3)}
    
    def load(self, dataID:str):
        """Loads the data from the specified stellar system.
        Raises Exception if the dataID is invalid.

        ----
        Arguments:
            dataID (str) -- The ID of the data to load.
        """
        self.dataID = dataID.upper()

        self.radius, self.mass = self.stellarMassRadiusDict[self.dataID] if self.dataID in self.stellarMassRadiusDict else (None, None)
        try:
            directory = self.systemsDirectoryDict[self.dataID]
            if directory.startswith("Data"):
                self.times, self.flux = np.loadtxt(directory, unpack=True, skiprows=3, usecols=[1,2])
            elif directory.startswith("NewData"):
                self.times, self.flux = np.loadtxt(directory, unpack=True, skiprows=136, usecols=[0,7])
                self.__removeNANFluxValues()

                lines = open(directory, 'r').readlines()
                #Parsing the /RADIUS file attribute
                self.radius = float(lines[44][10:-1])
                #Parsing and converting the /LLOG (stellar surface gravity) file attribute from log(cms^-2) to ms^-2
                surface_gravity = 10**(float(lines[40][8:-1]) - 2)
                self.mass = formulas.stellarMass(self.radius, surface_gravity)
        except KeyError:
            raise Exception("INVALID DATA ID: System not found.")

    def __removeNANFluxValues(self):
        nanValuesIndexes = [i for i, flux in enumerate(self.flux) if np.isnan(flux)]
        if len(nanValuesIndexes) > 0:
            self.times = np.delete(self.times, nanValuesIndexes)
            self.flux = np.delete(self.flux, nanValuesIndexes)

    def getData(self):
        return self.times, self.flux
    
    def getRadius(self):
        return self.radius

    def getMass(self):
        return self.mass
    
    def __iter__(self):
        for dataID in self.systemsDirectoryDict.keys():
            self.load(dataID)
            yield self