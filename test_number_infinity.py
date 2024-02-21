'''
Please ignore this file :)
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from data_analyser import DataAnalyser

analyser = DataAnalyser('kplr002853093')
analyser.plot('p')
plt.show()