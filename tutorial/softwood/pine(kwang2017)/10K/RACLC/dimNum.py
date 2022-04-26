import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
# import seaborn as sns
import os
import re
import csv
import math

def calcN2Lambda(T):
    lambdaN2 = 12.8 +0.05224 *T - 0.6482*10**(-6)*T**2 - 0.2765*10**(-9)*T**3
    lambdaN2 = lambdaN2/1000
    return lambdaN2

def convcCoeff(Nu_,lambdaN2_,dp_):
    ht_ = Nu_*lambdaN2_/dp_
    return ht_

def tConvc(rhop, Cpp, dp, ht):
    t = rhop * Cpp * dp/(6*ht)
    return t

def tCondc(rhop, Cpp, dp, lambdap):
    t = rhop * Cpp * dp**2/(36*lambdap)
    return t

def tRad(rhop, Cpp, dp, Tg, Tp, omega, sigma):
    t = rhop * Cpp * dp / ((6*omega*sigma)*(Tg+Tp)*(Tg**2+Tp**2))
    return t

def kPyroFunc(kCoeff, T, Ea):
    R = 8.3145 #[J mol^-1 K^-1]
    k = kCoeff * math.exp(-Ea/(R*T))
    return k

def kPyroFunc2(kCoeff, T, Ea):
    RR = 8314.5 #[J kmol^-1 K^-1]
    k = kCoeff * math.exp(-Ea*4184/(RR*T))
    return k

def NuCalc(ReD, Pr):
    Nu = 2 + 0.6*ReD**(1/2) * Pr**(1/3)
    return Nu

def NuCalc2(ReD, Pr):
    Nu = 2 + (0.4*ReD**(1/2)+0.06*ReD**(2/3))*Pr**(0.4)
    return Nu

def ReDCalc(U,L,nu):
    ReD = U*L/nu
    return ReD

def PrCalc(Cpg,mug,lambdag):
    Pr = Cpg*mug/lambdag
    return Pr

def muN2(T):
    mu0 = 1.663*10**(-5); T0 = 273; Smu = 107
    mu = mu0*(T/T0)**(3/2)*(T0+Smu)/(T+Smu)
    return mu

def rhoGAtm(MM,T):
    MMkg = MM/1000 # molar mass in kg/mol
    rhog = 10**5*MMkg/(8.314*T)
    return rhog

def CpN2(Tg):
    if 100<Tg<500:
        A = 28.98641; B = 1.853978; C = -9.647459; D = 16.63537; E = 0.000117;
    if 500<=Tg<2000:
        A = 19.50583; B = 19.88705; C = -8.598535; D = 1.369784; E = 0.527601;
    t = Tg/1000;
    CpN2 = A + B*t + C*t**2 + D*t**3 + E/(t**2) # J/mol*K
    CpN2 = CpN2*1000/28 # J/kg*K
    return CpN2


def dpSizeChange():
    Nu = 2; dp = 0.4*10**(-3); rhop = 300; Cpp = 1522;
    omega = 0.9; sigma = 5.670374 * 10**(-8)
    Tg = 1273; Tp = 300; lambdaN2 = calcN2Lambda(Tg); lambdap = 0.11
    kCoeff = 2.64*10**5; Ea = 105000
    kPyro = kPyroFunc(kCoeff, Tg, Ea)
    tPyro = 1/kPyro
    valEnum = 100 # Number of values

    # ht = convcCoeff(Nu,lambdaN2,dp)
    # tConvcExt = tConvc(rhop, Cpp, dp, ht)

    dpList = [(0.2*i+1)*dp for i in range(valEnum)]


    htList = [convcCoeff(Nu,lambdaN2,dpList[i]) for i in range(valEnum)]
    # print(htList)
    tConvcExtList = \
    [tConvc(rhop, Cpp, dpList[i], htList[i]) for i in range(valEnum)]
    tCondcInterList = \
    [tCondc(rhop, Cpp, dpList[i], lambdap) for i in range(valEnum)]
    tRadExtList = \
    [tRad(rhop, Cpp, dpList[i], Tg, Tp, omega, sigma) for i in range(valEnum)]
    # tRadExt = tRad(rhop, Cpp, dp, Tg, Tp, omega, sigma)
    BiList = \
    [tCondcInterList[i]/tConvcExtList[i] for i in range(valEnum)]


    tRadBytConv = [tRadExtList[i]/tConvcExtList[i] for i in range(valEnum)]
    tConvcBytPyro = [tConvcExtList[i]/tPyro for i in range(valEnum)]
    tCondBytPyro = [tCondcInterList[i]/tPyro for i in range(valEnum)]
    plt.plot(dpList,tRadBytConv,'r--',label="tRadBytConv")
    plt.plot(dpList,tConvcBytPyro,'b--',label=" tConvcBytPyro ")
    plt.plot(dpList,tCondBytPyro,'b--',label="tConvcBytPyro")
    plt.legend()
    # plt.show()
    plt.cla()
    # plt.plot(dpList,BiList,'k--',label="Bi")
    # print(BiList)

def readEffectiveK():
    root = "./0/"
    fileName = "weightedK_"
    data = pd.read_csv(os.path.join(root, fileName),
    header=None, names=['value']).values.tolist()
    data_split = [item[0].split(' ') for item in data]
    data_clean = [list(filter(None, item)) for item in data_split]
    T = [float(i[0]) for i in data_clean]
    weightedK_ = [float(i[1]) for i in data_clean]
    return T,weightedK_


def TChange():
    Xplot = [-1e6,1e6]
    pyLow = [0.1,0.1]
    pyHigh = [10,10]
    dpList = [0.177*10**(-3)]
    for dp in dpList:
        TCalc,weightedK_Calc = readEffectiveK()
        Nu = 2; rhop = 570; Cpp = 1522; Urel = 0.0148 #dp = 5*10**(-3);
        omega = 0.9; sigma = 5.670374 * 10**(-8)
        Tg = 300; Tp = 300; lambdap = 0.12
        dT = 10;
        # lambdaN2 = calcN2Lambda(Tg)
        valEnum = len(TCalc)#100 # Number of values
        TgList = TCalc#[(dT*i+Tg) for i in range(valEnum)]
        CpList = [3.867*T+103 for T in TgList]
        # print(CpList)
        lambdaN2List = [calcN2Lambda(value) for value in TgList]

        kCoeff = 2.64*10**5; Ea = 105000
        kPyroList = [kPyroFunc(kCoeff, TgList[i], Ea) for i in range(valEnum)]
        kHCE = 1.2500e11; EaHCE = 31400.00;
        KLIGH = 6.7000e12; EaLIGH = 37500.00;
        kPyroListHCE = [kPyroFunc2(kHCE, TgList[i], EaHCE) for i in range(valEnum)]
        kPyroListLIGH = [kPyroFunc2(KLIGH, TgList[i], EaLIGH) for i in range(valEnum)]
        # print(kPyroFunc2(kHCE,566,EaHCE)*0.004997882103)
        # print(weightedK_Calc)
        # tPyroList = [1/value for value in kPyroListHCE]
        # tPyroList = [1/value for value in kPyroList]
        # ht = convcCoeff(Nu,lambdaN2,dp)
        # tConvcExt = tConvc(rhop, Cpp, dp, ht)


        #calculate property of N2
        muN2List = [muN2(value) for value in TgList]
        rhoN2List = [rhoGAtm(28,value) for value in TgList]
        nuList = [muN2List[i]/rhoN2List[i] for i in range(valEnum)]
        CpN2List = [CpN2(value) for value in TgList]
        # calculate dimensionless number
        ReDList = [ReDCalc(Urel,dp,value) for value in nuList]
        PrList = \
        [PrCalc(CpN2List[i],muN2List[i],lambdaN2List[i]) for i in range(valEnum)]
        NuList = [NuCalc(ReDList[i],PrList[i]) for i in range(valEnum)]
        NuList2 = [NuCalc2(ReDList[i],PrList[i]) for i in range(valEnum)]
        NuDiff = [NuList[i]-NuList2[i] for i in range(valEnum)]
        # Nu2List =

        tPyroList = [1/value for value in weightedK_Calc]
        htList = [convcCoeff(NuList2[i],lambdaN2List[i],dp) for i in range(valEnum)]
        tConvcExtList = \
        [tConvc(rhop, CpList[i], dp, htList[i]) for i in range(valEnum)]
        tCondcInterList = \
        [tCondc(rhop, CpList[i], dp, lambdap) for i in range(valEnum)]
        tRadExtList = \
        [tRad(rhop, CpList[i], dp, TgList[i], Tp, omega, sigma) for i in range(valEnum)]
        tRadBytConv = [tRadExtList[i]/tConvcExtList[i] for i in range(valEnum)]
        tConvcBytPyro = [tPyroList[i]/tConvcExtList[i] for i in range(valEnum)]
        tCondBytPyro = [tPyroList[i]/tCondcInterList[i] for i in range(valEnum)]
        # BiList = [NuList[i]*lambdaN2List[i]/lambdap for i in range(valEnum)]
        BiList = [tCondBytPyro[i]/tConvcBytPyro[i] for i in range(valEnum)]
        # htList = [NuList[i]*lambdaN2List[i]/dp for i in range(valEnum)]
        # ht2List = [5.69+0.0098*T for T in TgList]
        # print(htList)
        # print(ht2List)
        # plt.plot(TgList,tRadBytConv,'r--',label="tRadBytConv")
        # print(TgList)

        plt.plot(TgList,tConvcBytPyro,'-',label="Py2"+" "+str(dp*1000)+"mm")
        plt.plot(TgList,tCondBytPyro,'--',label="Py1"+" "+str(dp*1000)+"mm")
        plt.plot(TgList,BiList,'r',label="Bi")
        plt.xlabel('Temperature (K)')
        plt.ylabel('Pyrolysis number')
        plt.yscale('log')
        plt.legend()
    plt.plot(Xplot,pyLow,'--')
    plt.plot(Xplot,pyHigh,'--')
    # plt.xlim([300, max(TCalc)])
    plt.xlim([min(TCalc), max(TCalc)])
    plt.savefig('PyByT.svg', format='svg',dpi=1200)
    plt.show()
    plt.cla()
    print('py1min: ', min(tConvcBytPyro))
    print('py2min: ', min(tCondBytPyro))

def biotNumDp():
    Urel = 0.0148; Tg = 1500; dp = 1*10**(-6);lambdap = 0.12
    valEnum = 2000 # Number of values
    dpList = [dp*(i+1) for i in range(valEnum)]

    muN2_ = muN2(Tg)
    rhoN2_ = rhoGAtm(28,Tg)
    nu = muN2_/rhoN2_
    CpN2_ = CpN2(Tg)
    lambdaN2_ = calcN2Lambda(Tg)
    ReDList = [ReDCalc(Urel,value,nu) for value in dpList]
    Pr_ = PrCalc(CpN2_,muN2_,lambdaN2_)
    NuList = [NuCalc(ReDList[i],Pr_) for i in range(valEnum)]
    BiList = [NuList[i]*lambdaN2_/lambdap for i in range(valEnum)]
    plt.plot(dpList,BiList,'r--',label="Bi")
    # plt.yscale('log')
    plt.legend()
    plt.show()
    plt.cla()


if __name__ == "__main__":
    # here we use constant gas temperature but increase the particle diameter
    # dpSizeChange()
    # # Here we are looking at the temperature influence to the
    # # charateristic values
    TChange()
    # biotNumDp()
