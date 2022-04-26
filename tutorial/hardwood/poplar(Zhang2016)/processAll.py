import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
# import seaborn as sns
import os
import re
import csv

def readCsvToList(fileName):
    with open(fileName) as tgaSim:
        tgaSimReader = csv.reader(tgaSim)
        T = list(); massLoss = list()
        for row in tgaSimReader:
            T.append(float(row[0]))
            massLoss.append(float(row[1]))
    return T,massLoss

T5K, massLoss5K = readCsvToList('5K/CRECK2017.csv')
T20K, massLoss20K = readCsvToList('20K/CRECK2017.csv')

plt.plot(T5K,massLoss5K,label='5K')
plt.plot(T20K,massLoss20K,label='20K')

plt.legend()
plt.show()
