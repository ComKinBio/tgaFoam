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


root = '../'
for folderName in os.listdir(root):
    if folderName=='postProcess':
        continue
    path=root+folderName
    filePath=path+'/'+folderName+'.csv'
    T, massLoss = readCsvToList(filePath)
    plt.plot(T,massLoss,label=folderName)
    # data = pd.read_csv(filePath,header=None)
    # data.plot(x =0, y=1)

plt.xlabel('Temperature (K)', fontweight ='bold', fontsize = 12)
plt.ylabel('mass loss (-)', fontweight ='bold', fontsize = 12)
# Read experimental data
expT, expMassLoss = readCsvToList('30K.csv')
plt.plot(expT,expMassLoss,'o',label='experiment')
plt.legend()
plt.xlim([min(expT),max(expT)])

plt.savefig('Lee2018Oak'+'.eps', bbox_inches='tight', dpi=400)
plt.show()
