import numpy as np
import os
import natsort
from datetime import date
import matplotlib.pyplot as plt

class System:
    def __init__(self, dataDir=None):
        self.fDataDir = dataDir
        self.fCodeData = dict()
        self.fSupplData = dict()

        self.fData = dict()
        self.fData["CODE"] = self.fCodeData
        self.fData["Suppl"] = self.fSupplData

    def __str__(self):
        return self.fData

    # --- * --- Set methods
    def setDataDir(self, dataDir):
        self.fDataDir = dataDir

    def addCodeData(self, *args):
        for codeData in args:
            codeName = codeData.getCodeName()
            self.fCodeData[codeName] = codeData
        
    def addSupplementalData(self, *args):
        for supplementalData in args:
            supplementalDataName = supplementalData.getCodeName()
            self.fSupplData[supplementalDataName] = supplementalData

    # set position List 
    def setPositionList(self, positionList):
        for data in self.fCodeData.items():
            data = data[1]
            data.setPositionList(positionList)

    def setPostureList(self, postureList):
        for data in self.fCodeData.items():
            data = data[1]
            data.setPostureList(postureList)
    
    # --- * --- Get methods
    def getDataDir(self):
        return self.fDataDir

    def getCodeData(self):
        return self.fCodeData

    def getSupplData(self):
        return self.fSupplData

    def getRawData(self):
        Data = dict()
        for data in self.fCodeData.items():
            codeName = data[0]
            data = data[1]
            Data[codeName] = data.getRawData()
        return Data

    def getEffDoseData(self):
        Data = dict()
        for data in self.fCodeData.items():
            codeName = data[0]
            data = data[1]
            Data[codeName] = data.getEffDose()
        return Data

    def getEDRsData(self):
        Data = dict()
        for data in self.fCodeData.items():
            codeName = data[0]
            data = data[1]
            Data[codeName] = data.getEDRs()
        return Data

    def getOrganContr(self):
        Data = dict()
        for data in self.fCodeData.items():
            codeName = data[0]
            data = data[1]
            Data[codeName] = data.getOrganContr()
        return Data

    def getAnnualDose(self):
        Data = dict()
        for data in self.fCodeData.items():
            codeName = data[0]
            data = data[1]
            Data[codeName] = data.getAnnualDose()
        return Data

    # --- * --- Computing methods
    def calcEffectiveDose(self):
        for data in self.fCodeData.items():
            data = data[1]
            suppleData = self.getSupplData()["TWF"]
            data.calcEffDose(suppleData.getData())

    def writeEffDose(self):
        for data in self.fCodeData.items():
            data = data[1]
            data.writeResultEffDose()

    def writeOrganEffDose(self):
        for data in self.fCodeData.items():
            data = data[1]
            data.writeOrganEffDose()

    def computeEDR(self):
        for data in self.fCodeData.items():
            data = data[1]
            data.computeEDR()

    def computeOrganEDR(self):
        for data in self.fCodeData.items():
            data = data[1]
            data.convertOrganEffDoseContr()
    
    # --- * --- Plotting Area
    def init_plot(self, imageFolderPath: str):
        for data in self.fCodeData.items():
            data = data[1]
            data.init_plot(imageFolderPath=imageFolderPath)
        self.fImageFolderPath = imageFolderPath

    def drawBasicFigures(self):
        for data in self.fCodeData.items():
            data = data[1]
            data.drawBasicFigures()

    def drawCombinedFigures(self):
        today = str(date.today())
        imageFolderDir = os.path.join(self.fImageFolderPath, today)
        saveFolderPath = os.path.join(imageFolderDir, "NoTrend")

        def reset(self):
            figure = plt.figure(figsize=(12, 8))
            self.figure = figure
            self.fmt = ["o", "^", "s"]
            ax = figure.add_subplot(111)
            self.ax = ax

        reset(self)
        def init(self):
            self.fOptions = {"legendFontSize": 20,
                            "xLabeltxt": "Position (cm)",
                            "yLabeltxt": "Effective Dose Rate (\u03BCSv/s)",
                            "labelFontSize": 15,
                            "xticksFontSize": 15,
                            "xScale": "linear",
                            "yScale": "log",
                            "title": "Effective Dose Rate for Operating Worker",
                            "titleFontSize": 20,
                            "grid": {"alpha": 0.7, "color": "gray", "linestyle": "--"},
                            "xLim": [-160, 110],
                            "yLim": [5.e-07, 1]
                            }
        init(self)
        
        def applyOption(self):
            self.ax.legend(fontsize = self.fOptions["legendFontSize"])
            self.ax.set_xlabel(self.fOptions["xLabeltxt"], fontsize = self.fOptions["labelFontSize"])
            self.ax.set_ylabel(self.fOptions["yLabeltxt"], fontsize = self.fOptions["labelFontSize"])
            plt.xticks(fontsize = self.fOptions["xticksFontSize"])
            plt.yticks(fontsize = self.fOptions["xticksFontSize"])
            self.ax.set_xscale(self.fOptions["xScale"])
            self.ax.set_yscale(self.fOptions["yScale"])
            self.ax.set_title(self.fOptions["title"], fontsize = self.fOptions["titleFontSize"])
            self.ax.grid(alpha = self.fOptions["grid"]["alpha"], 
                    color = self.fOptions["grid"]["color"], 
                    linestyle = self.fOptions["grid"]["linestyle"]
                    )
            self.ax.set_xlim(self.fOptions["xLim"][0], self.fOptions["xLim"][1])
            self.ax.set_ylim(self.fOptions["yLim"][0], self.fOptions["yLim"][1])
        
        keys = list(self.fCodeData.keys())
        postureList = self.fCodeData[keys[0]].getPostureList()
        for posture in postureList:
            self.fOptions["title"] = f"Effective Dose Rate for {posture} Worker"
    
            for key in keys:
                data = self.fCodeData[key]
                (x,y,yerr,codeName) = data.exportPlotting(posture)
                if np.isnan(y)[0] == True:
                    continue
                self.ax.errorbar(x = x, y = y, yerr = yerr, label = codeName, fmt = "o")
            
            applyOption(self)
            savePath = os.path.join(saveFolderPath, f"{posture}_LOG.svg")
            if not os.path.isdir(saveFolderPath):
                    os.makedirs(saveFolderPath)

            plt.savefig(savePath, dpi=300)
            plt.close()
            reset(self)


    def exportCoef(self, fittool):
        for posture in fittool.fPostureList:
            for data in self.fCodeData.items():
                data = data[1]
                data.getCoef(posture=posture, data = fittool.exportCoef(posture=posture, codeName=data.getCodeName()))
