import matplotlib
import math
import csv
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

def searchPhaseNames(phases_, phaseName_, splitElement_):
    nameString_ = phases_.split(phaseName_[splitElement_])[1]\
               .split(phaseName_[splitElement_+1])[0].split('{')[1] .replace("}", "")
    namesList_ = list()
    names = nameString_.split(';')
    for i in names:
        namesList_.append(i.split())

    namesList_ = list(filter(None, namesList_))
    return namesList_

def makeNameList(phaseNameList_, state):
    outputList = list()
    for name in phaseNameList_:
        outputList.append("Y" + name[0] + state)
    return outputList


time = 300
#fileNames = list()
fileNamesList = list()
fileNamesFile = open("constant/coalCloud1Properties")

# read all contents from the file above
contentCoalClouds = fileNamesFile.read()
phases = contentCoalClouds.split('phases')[1]

phasesName = ['gas', 'liquid', 'solid', '}'] # the { is used to separate

for splitElement in range(len(phasesName)-1):
    globals()[phasesName[splitElement]] \
        = searchPhaseNames(phases, phasesName, splitElement)

gasFileNameList = makeNameList(gas, "(g)")
solidFileNameList = makeNameList(solid, "(s)")
liquidFileNameList = makeNameList(liquid, "(l)")

# Read file in this time interval
endTime = 21000
timeRange = range(30,endTime,30)

# Read temperature profile
T = list()
for time in timeRange:
    fileName = str(time) + "/lagrangian/coalCloud1/T"
    TOpen = open(fileName).readlines()[17]
    Tval = float(TOpen.split(')')[0].replace("1(", ""))
    T.append(Tval)

# plot the results for lagrangian solid phase
for solidFileNames in range(len(solidFileNameList)):
    variables = list()
    for time in timeRange:
        fileName = str(time) + "/lagrangian/coalCloud1/" +solidFileNameList[solidFileNames]
        globals()[solidFileNameList[solidFileNames]] = open(fileName)
        content = globals()[solidFileNameList[solidFileNames]].readlines()
        line = content[17]
        var = float(line.split(')')[0].replace("1(", ""))
        variables.append(var)
        globals()[solidFileNameList[solidFileNames]].close()
    plt.plot(T, variables, label=solidFileNameList[solidFileNames])

plt.xlabel('Temperature (K)')
plt.ylabel('Species fractions (-)')
plt.title('Species fractions')
plt.ylim([0, 1])
plt.xlim([300, 1200])
plt.legend(ncol=1, fontsize=7, bbox_to_anchor=(1, 1), loc='best', borderaxespad=1)
plt.savefig('speciesFraction.png', bbox_inches='tight', dpi=400)
plt.show()

# Read the mass0
mass0Open = open("30/lagrangian/coalCloud1/mass0")
mass0 = float(mass0Open.readlines()[17].split(')')[0].replace("1(", ""))

# The ash mass
ashOpen = open("30/lagrangian/coalCloud1/Yash(s)")
mAsh = float(ashOpen.readlines()[17].split(')')[0].replace("1(", "")) * mass0
print(mAsh)

mass0 = mass0 - mAsh

expT=list()
expMassLoss=list()
# Read experimental data from csv file
with open('2015data.csv') as expFile:
    expFileReader = csv.reader(expFile)
    for row in expFileReader:
        expT.append(float(row[0]))
        expMassLoss.append(float(row[1]))

# From celsius degree to K
expT = [expT_+273 for expT_ in expT]

# The expected experiment result/ previous numerical results
numT = [20.64516129, 96.77419355, 165.1612903, 237.4193548, 273.5483871,
        294.8387097, 318.0645161, 307.0967742, 325.1612903, 338.0645161,
        352.2580645, 385.8064516, 415.483871, 444.516129, 467.0967742,
        505.8064516, 552.2580645, 594.8387097, 638.7096774, 677.4193548, 700]


numT = [expT_+273 for expT_ in numT]

numMassLoss = [1, 0.996070727, 0.99410609, 0.956777996, 0.882121807,
               0.815324165, 0.648330059, 0.746561886, 0.565815324,
               0.418467583, 0.326129666, 0.286836935, 0.275049116,
               0.259332024, 0.249508841, 0.239685658, 0.239685658,
               0.233791749, 0.233791749, 0.229862475, 0.223968566]
               

# plot the total mass of the solid
massit = list()
for time in timeRange:
    rhoFileName = str(time) + "/lagrangian/coalCloud1/rho"
    dFileName = str(time) + "/lagrangian/coalCloud1/d"

    rhoOpen = open(rhoFileName)
    rhoContent = rhoOpen.readlines()[17]
    rhoVal = float(rhoContent.split(')')[0].replace("1(", ""))

    dOpen = open(dFileName)
    dContent = dOpen.readlines()[17]
    dVal = float(dContent.split(')')[0].replace("1(", ""))

    # the mass change w.r.t time
    massit.append(((1./6.)*math.pi*rhoVal*(dVal**3)- mAsh)/mass0)

# Compare the numerical results with the experiment
plt.plot(T, massit, label="numerical mass loss")
plt.plot(expT, expMassLoss, 'x', label="experimental mass loss")
plt.plot(numT, numMassLoss, label="old simulation mass loss")
plt.xlabel('Temperature (K)')
plt.ylabel('Mass loss (-)')
plt.title('TGA Curve')
plt.legend(ncol=1, fontsize=7)
plt.legend(bbox_to_anchor=(1, 1), loc='best', borderaxespad=0.5)
plt.savefig('massFraction.png', dpi=400)
plt.show()


