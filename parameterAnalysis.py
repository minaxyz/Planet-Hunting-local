from data_analyser import DataAnalyser
from data_handler import LocalDataHandler
for handler in LocalDataHandler():
    print(handler.getRadius())