# Planet Hunting with Python

## Introduction

This is a repository that builds on provided code from QMUL. We analyse data sent from NASA's Kepler spacecraft 
## Next steps

- plot models to data, optimise accuracy of the models
- draw conclusions about them, analyse how big they are, how many planets in a system perhaps

## What this repository contains

- Module for finding planet parameters
- Module for plotting piece wise function (line of best fit) | in progress
- Module for Chi-squared test
- Module for extension tasks | planned

## List of files

## utils.py
- Imports
floor and ceil from math
- Functions

binarySearch()

    * Parameters : array, element
    * finds an item in a sorted array using a binary search method
    * O(log n) time complexity
    * returns the location of the element in the array

interpolationSearch()

    * Parameters : array, element
    * finds an item in a sorted array using a interpolation search method
    * O(log log n) time complexity - an improvement on binary search
    * returns the location of the element in the array

## utils_test.py
- Imports
the information / functions from utils.py
perf_counter from time
random
numpy
- Functions

benchmarkArrayFunc()

    * Parameters : size (size of the array created), n (number of times each search is performed), funcs(the functions to be tested)
    * compares the average time of different inputted functions
    * O(n²)
    * returns the name of the function and average time in milliseconds for each of the functions inputted

## test.py
- Imports
DataAnalyser from data_analyser
perf_counter from time
- Classes
   - Timing()

- Functions (in Timing())

__init__()

   * Parameters : sinceLastOut(set to be true on defult), ms(set to be true on defult)
   * self.s (the start time of when the class is instantiated)
   * self.last (set to equal self.s)
   * self.sinceLastOut (set to equal sinceLastOut)
   * self.unit (value will be 1000 if ms is true, else the value will be 1)
   * self.unitPrefix (value will be 'ms' if ms is true, else the value will be 's')

out()

   * Parameters : self, label
   * prints the value of label then the time elapsed since the class was instantiated in the unit chosen and with the value in unitPrefix concatinated onto the end
   * if sinceLastOut is still true then self.last will be given the value of the current time

totalOut()

   * Parameters : self, label(set to be an empty string on defalt)
   * if self.sinceLastOut has the value False then self.out() will be called with the parameter label
   * then, regardless of the outcome, it prints 'Total:' followed by the time elapsed since the class was instantiated in the unit chosen and with the value in unitPrefix concatinated onto the end

- Functions

timedTest()

   * Parameters : dataID, plotType(set to be an empty string on defalt)
   * prints the dataID followed by 'results:'
   * the the class Timing is instantiatied under the variable t and the class DataAnalyser, from ##data_analyser is instantiatied under the variable analyser
   * the program then outputs the time taken for initialisation and then calls the function getOrbitalPeriod() in analyser and outputs the time taken for that function to perform
   * after, the prorgam calls the function getModel() from analyser and outputs the values of m.min and m.max. then the program outputs the time taken for that function to perform
   * then the program prints the value of period and calls the function totalOut() in t. After the data is plotted from the fuction plot() from analyser

## pwlf.py
piecewise linear library -- not our code but referred to in other parts of the project

## data_analyser.py
