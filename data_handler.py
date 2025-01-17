import formulas
from abc import abstractmethod
from PyPDF2 import PdfReader
import numpy as np
import os

LOCAL_DATA_FOLDERS = ["Data", "NewData"]

def setDirectories(directories):
    """Specifies the directories in which to search for lightcurve files (any .tbl or .dat file).
    By default "Data" and "NewData".

    Arguments:
        directories (iterable) - Directories to search for valid .dat or .tbl files (eg. ["Lightcurves", "Lightcurves/Kepler"]). 
    """
    global LOCAL_DATA_FOLDERS
    LOCAL_DATA_FOLDERS = directories

class AbstractDataHandler():
    """Abstract data handler class for the creation of data handlers which are compatible with data analyser.
    All data handlers must support `getData`, `getRadius`, and `getMass` methods, and optionally the `__iter__` method to allow for
    data analyser to become iterable as well.
    """
    def __init__(self, dataID):
        """Loads the data from the specified stellar system.

        ----------
        Arguments:
            dataID (str) - The ID of the data to load.
        """
        self.times, self.flux, self.radius, self.mass = None, None, None, None 

    def getData(self):
        """Returns the flux against time data of the stellar system.
        
        ----------
        Returns:
            tuple (times, flux):
                times (np.array) - Contains the times for the corresponding flux recordings -- Units: `Days`.

                flux (np.array) - Contains the recordings of the intensity of the star's light. Data must be normalised
                (i.e. flux must be a fraction of the mean).

            *Note: The arrays returned do not contain NaN values.
        """
        return self.times, self.flux
    
    def getRadius(self):
        """Returns the fetching of the star's radius.
        
        ----------
        Returns:
            radius (float) - The radius of the star -- Units: `solar radii`.
        """
        return self.radius
    
    def getMass(self):
        """For the fetching of the star's mass.

        ----------
        Returns:
            mass (float) - The mass of the star -- Units: `solar masses`.
        """
        return self.mass

class LocalDataHandler(AbstractDataHandler):
    """Built-in fetching of data from specified local files which follow the old or new nasa exoplanet archive data format.
    (By default it is assumed that files follow the new format, so if a file follows the old format, 
    the file's name needs to start with KIC)
    LocalDataHandler is iterable.
    """
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
        self.systemsDirectoryDict = {}
        for folderName in LOCAL_DATA_FOLDERS: 
            self.systemsDirectoryDict |= {
                file.split('.')[0].split('_')[0].upper():folderName + '/' + file 
                for file in os.listdir(folderName) 
                if file.endswith((".tbl",".dat")) 
                and "phaseFold" not in file and "TIC" not in file}
    
    def __initialiseStellarMassAndRadius(self):
        self.stellarMassRadiusDict = None
        try:
            f = PdfReader(open("Data/Stellar_Mass_Radius.pdf", 'rb'))
            tableData = f.pages[0].extract_text().split()[9:]
            self.stellarMassRadiusDict = {tableData[i]:(float(tableData[i+2]),float(tableData[i+1]))  for i in range(0,len(tableData),3)}
        except Exception:
            pass
    
    def load(self, dataID:str):
        """Loads the data from the specified stellar system.
        Raises Exception if the dataID is invalid.

        ----------
        Arguments:
            dataID (str) -- The ID of the data to load.
        """
        self.dataID = dataID.upper()

        self.radius, self.mass = self.stellarMassRadiusDict[self.dataID] if self.dataID in self.stellarMassRadiusDict else (None, None)
        try:
            directory = self.systemsDirectoryDict[self.dataID]
            if "KIC" in directory.upper():
                self.times, self.flux = np.loadtxt(directory, unpack=True, skiprows=3, usecols=[1,2])
            else:
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

    def __iter__(self):
        for dataID in self.systemsDirectoryDict.keys():
            self.load(dataID)
            yield self