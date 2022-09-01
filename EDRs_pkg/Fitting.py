import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import natsort
import os
from scipy.optimize import curve_fit
from datetime import date

class Fitting:
    def __init__(self):
        # data containers
        self.fData = dict()

        self.fEDRs = dict()
        self.fRavelEDRs = dict()

        self.fOrganContr = dict()

        # interpolation
        self.fCoef = dict()

        # supplementary data
        self.fPositionList = []
        self.fPostureList = []
        self.fCodeNames = []

    def addData(self, *args):
        for dataClass in args:
            dataName = dataClass.getCodeName()
            self.fCodeNames.append(dataName)
            self.fPositionList = dataClass.getPositionList()
            self.fPostureList = dataClass.getPostureList()
            self.fData[dataName] = dataClass

        for posture in self.fPostureList:
            ### EDRs / Organ Contr.
            self.fEDRs[posture] = dict()
            self.fOrganContr[posture] = dict()
            self.fCoef[posture] = dict()
            for codeName in self.fCodeNames:
                self.fEDRs[posture][codeName] = pd.DataFrame(self.fData[codeName].getDF_EDRs(posture))
                self.fOrganContr[posture][codeName] = pd.DataFrame(self.fData[codeName].getDF_OrganContr(posture))
                self.fCoef[posture][codeName] = dict()
        
        self.ravelDF()
    
    def ravelDF(self):
        for posture in self.fPostureList:
            self.fRavelEDRs[posture] = pd.DataFrame()
            for codeName in self.fCodeNames:
                self.fEDRs[posture][codeName].rename(columns= {"EDR": f"{codeName}_EDR",
                                                            "EDR_RelErr": f"{codeName}_RelErr"}, inplace=True)
                self.fRavelEDRs[posture] = pd.concat([self.fRavelEDRs[posture], self.fEDRs[posture][codeName]], axis=1)
    
    # --- interpolation function
    def getCoef(self, x, y, yerr, field: str):
        """ interpolate points for under ground, field: ground, under, and total"""
        
        def func(x, a, b):
            return a + b * x
        
        field = field.lower()

        if field == "ground":
            newx = x[:10]
            newy = y[:10]
            newyerr = yerr[:10]
            
        elif field == "under":
            newx = x[9:]
            newy = y[9:]
            newyerr = yerr[9:]
            
        elif field == "total":
            newx = x
            newy = y
            newyerr = yerr
            
        logy = np.log(newy)
        
        popt, pcov = curve_fit(func, newx, logy, method="lm")
        
        # get r squared
        c = np.exp(popt[0]) * np.exp(popt[1] * newx)
        correlation_matrix = np.corrcoef(newy, c)
        correlation_xy = correlation_matrix[0, 1]
        r_squared = correlation_xy**2
        
        print(np.exp(popt), r_squared)

        return (newx, newy, newyerr, (np.exp(popt[0]), popt[1]), r_squared)

    def getCoefCombined(self, x, y, yerr, field: str, sZPos = 98.92301, sYPos = 44.92301, pdap = False): 
        """ get coef. for ground state, """
        def func(x, a, b, c, d):
            # 44.92301 , 98.92301
            r = np.sqrt(sYPos ** 2 + (sZPos - x) ** 2)
            return a * x**3 + b * x**2 + c * x**1 + d

        # Ground
        if field == "ground":
            if pdap:
                newx = x[5:10]
                newy = y[5:10]
                newyerr = yerr[5:10]
            else:
                newx = x[:10]
                newy = y[:10]
                newyerr = yerr[:10]
            
        elif field == "under":
            newx = x[10:]
            newy = y[10:]
            newyerr = yerr[10:]
            
        elif field == "total":
            newx = x
            newy = y
            newyerr = yerr
            
        p0 = (100, -1.0, 50, 0.1)
        popt, pcov = curve_fit(func, newx, newy, p0=p0, method="lm", maxfev=100000)
        
        # get r squared
        c = func(newx, *popt)
        correlation_matrix = np.corrcoef(newy, c)
        correlation_xy = correlation_matrix[0, 1]
        r_squared = correlation_xy**2
        
        print(*popt, "\t R^2: ",r_squared)
        return (newx, newy, newyerr, popt, r_squared)

    def interpolate(self):
        sYpos = 44.92301
        sZpos = 98.92301
        pdap = False

        # Ground
        for posture in self.fPostureList:

            if posture == "sp":
                sYpos = 124.92301

            elif posture == "sq":
                sYpos = 48.7
                sZpos = 66.3
                pdap = True

            elif posture == "bd":
                sYpos = 42.4
                sZpos = 90.3
                pdap = True

            for codeName in self.fCodeNames:
                x = np.array(self.fPositionList)
                y = np.array(self.fRavelEDRs[posture][f"{codeName}_EDR"])
                if np.isnan(y)[0] == True:
                    continue
                yerr = np.array(self.fRavelEDRs[posture][f"{codeName}_RelErr"]) * y
                print(posture, codeName)
                self.fCoef[posture][codeName]["ground"] = self.getCoefCombined(x = x, y = y, yerr = yerr, \
                    field = "ground", sYPos=sYpos, sZPos=sZpos, pdap=pdap)

        # Under
        for posture in self.fPostureList:
            for codeName in self.fCodeNames:
                x = np.array(self.fPositionList)
                y = np.array(self.fRavelEDRs[posture][f"{codeName}_EDR"])
                if np.isnan(y)[0] == True:
                    continue
                yerr = np.array(self.fRavelEDRs[posture][f"{codeName}_RelErr"]) * y
                print(posture, codeName)
                self.fCoef[posture][codeName]["under"] = self.getCoef(x = x, y = y, yerr = yerr, field = "under")

    def exportCoef(self, posture, codeName):
        return self.fCoef[posture][codeName]