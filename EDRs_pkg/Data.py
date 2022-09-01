import os
import natsort
import pandas as pd
import numpy as np
import math
import pickle
from datetime import date
import matplotlib.pyplot as plt

class Data:
    def __init__(self, filePath: str = "", codeType: str = ""):
        # --- * --- Data information
        self.fDirPath = filePath
        self.fFilePath = filePath
        self.fTypeOfCode_lower = codeType.lower()
        self.fTypeOfCode = codeType.upper()

        # --- * --- Data Variables
        self.fData = dict()

    def __str__(self):
        string = ""
        for key in self.fData.keys():
            string += f"{key}: {self.fData[key]}\n"
        return string
        
    # --- * --- Set methods
    def setFilePath(self, filePath):
        self.fFilePath = filePath

    def setCodeName(self, codeName: str):
        self.fTypeOfCode = codeName.upper()

    # --- * --- Get methods
    def getFilePath(self):
        return self.fFilePath
    
    def getCodeName(self):
        return self.fTypeOfCode
        
class MonteCarloData(Data):
    def __init__(self, *args, **kwargs):
        """ Inheritance of Partents' classes"""
        super().__init__(*args, **kwargs)
        # super(Plotting).__init__(*args, **kwargs)

        # --- * --- Data information
        self.fFilePath = os.path.join(self.fFilePath, self.fTypeOfCode_lower)
        self.fOrganIDs = []
        self.fPostureList = []
        self.fPositionList = []
    
        # --- * --- Data Variables
        self.fLatestFiles = list()
        self.fRawDataDic = dict()
        self.fEffectiveDose = dict()
        self.fEDRs = dict()
        self.fOrganContr = dict()
        self.initialize()
        self.fCoef = dict()

    # --- * --- Set methods
    def setLatestFiles(self, latestFiles: list):
        self.fLatestFiles = latestFiles

    def setData(self, data):
        self.fData = data

    def setPositionList(self, data: list):
        self.fPositionList = data

    def setPostureList(self, data: list):
        self.fPostureList = data

    def setDataFrame_EDRs(self):
        self.fDF_EDRs = pd.DataFrame(self.fEDRs, index=["EDR", "EDR_RelErr"])
        self.fDF_EDRs = self.fDF_EDRs.transpose()

    def setDataFram_OrganContr(self):
        self.fDF_OrganContr = pd.DataFrame(self.fOrganContr)
        self.fDF_OrganContr = self.fDF_OrganContr.transpose()

    # --- * --- Get methods
    def getLatestFiles(self):
        return self.fLatestFiles

    def getData(self):
        return self.fData

    def getRawData(self):
        return self.fRawDataDic
    
    def getEffDose(self):
        return self.fEffectiveDose

    def getEDRs(self):
        return self.fEDRs
    
    def getOrganContr(self):
        return self.fOrganContr

    def getAnnualDose(self):
        return self.fAnnualDose

    def getPositionList(self):
        return self.fPositionList

    def getPostureList(self):
        return self.fPostureList
    
    def getDF_EDRs(self, workType = ""):
        if workType != "":
            workType = workType.upper()
            return self.fDF_EDRs[workType]
        else:
            return self.fDF_EDRs

    def getDF_OrganContr(self, workType = ""):
        if workType != "":
            workType = workType.upper()
            return self.fDF_OrganContr[workType]
        else:
            return self.fDF_OrganContr
    
    # --- * --- Computing methods
    def initialize(self):
        # sub functions
        self.updateLatestFiles()

    def _getLatestFolderPath(self, filePath):
        LatestFolderPath = os.path.join(filePath, natsort.natsorted([folder for folder in 
                                        (os.listdir(filePath)) if folder.startswith("202")])[-1])
        # check existance of .DS_Store file in list
        for i in range(len(LatestFolderPath)):
            if ".DS_Store" in LatestFolderPath[i]:
                del LatestFolderPath[i]

        return LatestFolderPath
    
    def _getLatestFiles(self, latestFolderPath):
        LatestFiles = [os.path.join(latestFolderPath, _path)
                                    for _path in natsort.natsorted(os.listdir(latestFolderPath), reverse=True)]
        # check existance of .DS_Store file in list
        for fileName in LatestFiles:
            if ".DS_Store" in fileName:
                LatestFiles.remove(fileName)

        return LatestFiles
        
    def updateLatestFiles(self):
        rawFileDir = self.getFilePath()
        self.setLatestFiles(
                self._getLatestFiles(
                    self._getLatestFolderPath(rawFileDir) ) )
    
    def readData(self, isCombined = False):
        if self.fTypeOfCode == "GEANT4":
            self._readData()
        elif self.fTypeOfCode == "MCNP6":
            self._readData_MCNP6(isCombined=isCombined)
        elif self.fTypeOfCode == "PHITS":
            self._readData_PHITS(isCombined=isCombined)

    # Read fRawDataDic
    def _readData(self):
        for fileName in self.fLatestFiles:
            self.fRawDataDic[fileName] = dict()
            self.fRawDataDic[fileName] = dict()
            numberOfOrgans = 187
            with open(fileName) as f:
                line = f.readline()
                while True:
                    try:
                        line = f.readline().split()
                        if line[0] == "100":
                            for i in range(numberOfOrgans):
                                organID = line[0]
                                equivDose = float(line[-2])
                                absError = float(line[-1]) * equivDose
                                self.fRawDataDic[fileName][organID] = [equivDose, absError]
                                line = f.readline().split()
                        line = f.readline().split()
                    except:
                        break
        self.setData(self.fRawDataDic)
        
    def _readData_MCNP6(self, isCombined):
        def makeCombinedFiles():
            nps1e6 = "1000000"; nps1e7 = "10000000"; nps2e5 = "200000"; nps3e7 = "30000000"; nps5e5 = '500000'; nps25e5 = "250000"; nps5e6 ="5000000"; nps1e5 = "100000"
            npsList = [nps1e6, nps1e7, nps2e5, nps3e7, nps5e5, nps25e5, nps5e6, nps1e5]
            ResultDict = dict() # ResultDict[filename] =[(organID, equivDose, relativeError), ... (...)]
            # get data 

            for folderPath in self.fLatestFiles:
                
                ResultDict[folderPath] = dict()
                # equip..0 - equip..9
                lst = [os.path.join(folderPath, _path) for _path in (os.listdir(folderPath))]
                for i in range(len(lst)):
                    if ".DS_Store" in lst[i]:
                        print(i, len(lst))
                        del lst[i]
                        break
                        
                # print(folderPath)
                # print(lst)
                for num in range(len(lst)):
                    temp = dict() # temp[cellID] = [totalDose, relative Error]
                    # print (f"{num}th iteration started |n")
                    f = open(lst[num], "r") # open file with that name
                    # print (lst[num], "is opend")
                    while True:
                        # read lines and split it by space so that I can get the list of each of words in lines
                        line = f.readline();

                        if not line: break; # if the line reaches to eof, then break
                        list_line = line.split()

                        if len(list_line)==0: # skip the \n lines
                            continue    
                        # neutron Parts
                        # find the oragan ID and mass.  
                        if (list_line[0] == "1tally") and (len(list_line)==5 ) and \
                            (list_line[4] in npsList):
                            # ( list_line[4]==str(nps1e6) or list_line[4]==str(nps1e7)):

                            # find particle definition
                            while not line.strip().startswith("particle(s):"):
                                line = str(f.readline()) 
                            #set the ptl to particle name
                            ptl = (line.split())[1]

                            # find cell ID (organ ID)
                            while not line.strip().startswith("cell:"):
                                line = str(f.readline())
                            #set the cellID to organID
                            cellID = (line.split())[1]
                            
                            # check there a temp[cellID]
                            try: 
                                temp[cellID]
                            except:
                                temp[cellID] = [0., 0.] # temp[cellID] = [totalDose, relative Error]

                            # neutrons, dE and dF card used
                            if ptl == "neutrons":
                                while not line.strip().startswith("cell"):
                                    line = str(f.readline()) # cell:
                                line = str(f.readline()); line = str(f.readline()); line = str(f.readline()); line = str(f.readline()) # data line
                                equivDose = float(line.split()[0]); absDoseUncertainity = equivDose * float(line.split()[1])
                                temp[cellID][0] += equivDose; temp[cellID][1] = math.sqrt(temp[cellID][1]**2 + absDoseUncertainity**2)

                            # other particles, EM card used 
                            elif ptl == "photons" or ptl == "protons" or ptl == "alphas":
                                while not line.strip().startswith("energy"):
                                    line = str(f.readline())
                                line = str(f.readline()) # data line
                                equivDose = float(line.split()[1]); absDoseUncertainity = equivDose * float(line.split()[2])
                                temp[cellID][0] += equivDose; temp[cellID][1] = math.sqrt(temp[cellID][1]**2 + absDoseUncertainity**2)                        
                    f.close()
                    ## Set ResultDict[folderPath]
                    for organID in temp.keys():
                        totalEffDose = temp[organID][0]; totalEffDoseRelUncer = temp[organID][1]
                        # print(lst[num])
                        # print(organID, temp[organID])
                        try:
                            ResultDict[folderPath][organID]
                        except:
                            ResultDict[folderPath][organID] = [(totalEffDose, totalEffDoseRelUncer)]
                        else:
                            ResultDict[folderPath][organID].append( (totalEffDose, totalEffDoseRelUncer) )
                
                ## Error Propagation
                for organID in ResultDict[folderPath].keys():
                    effDoselst = [_effDose[0] for _effDose in ResultDict[folderPath][organID]]
                    effDoseAbsErrorlst = [_effDose[1] * _effDose[1] for _effDose in ResultDict[folderPath][organID]]
                    # print(organID, effDoselst, effDoseAbsErrorlst)
                    ResultDict[folderPath][organID][0] = sum(effDoselst) / len(lst)
                    # print(ResultDict[folderPath][organID])
                    if ResultDict[folderPath][organID][0] != 0:
                        try:
                            ResultDict[folderPath][organID][1] = (math.sqrt(sum(effDoseAbsErrorlst)) / len(lst))/ ResultDict[folderPath][organID][0]
                        except:
                            ResultDict[folderPath][organID].append((math.sqrt(sum(effDoseAbsErrorlst)) / len(lst))/ ResultDict[folderPath][organID][0])
                    else:
                        # print(lst[num])
                        # print(organID)
                        # print(ResultDict[folderPath][organID])
                        try:
                            ResultDict[folderPath][organID][1] = 0
                        except:
                            ResultDict[folderPath][organID].append(0)
                
                
                # print(ResultDict[folderPath])
                ## set resultFolder
                today = date.today()
                ## Data file Dir
                dataFileDir = self.fDirPath
                ## output files
                mcnp6FolderDir = os.path.join(dataFileDir, "mcnp6")
                ## resultFolderPath
                resultFolderPath = os.path.join(mcnp6FolderDir, f"combined_files/{today}")
                ## check dir 
                if not os.path.isdir(resultFolderPath):
                    os.makedirs(resultFolderPath)
                ## set the resultFileName
                # for files in mcnp6DataFilesList:

                dataType = folderPath.split("/")[-1]

                # for files in phitsDataFilesList:
                postureLst = ["op", "sp"]
                positionLst = [str(100 - pos*10) for pos in range(26)]
                positionLst = positionLst + ["7", "5", "3"]
                for posture in postureLst:
                    for position in positionLst:
                        text = f"{posture}_{position}"
                        if text == dataType:
                            resultFile = os.path.join(resultFolderPath, f"{posture}_{position}.out")
                
                g=open(resultFile, "w")
                g.write("organID \t equiv_Dose (Sv/source) \t relative error \n")
                for organID in ResultDict[folderPath].keys():
                    # '%-10s' % 'python'
                    sentence = '%-14s' % (str(organID)) + " | " + "{:.6e}".format(ResultDict[folderPath][organID][0])\
                        + "\t" +  "{:.6e}".format(ResultDict[folderPath][organID][1])
                    g.write(sentence + "\n")
                g.close()        
            self.fRawDataDic = ResultDict
            self.setData(self.fRawDataDic)

        if isCombined:
            makeCombinedFiles()
            
        self.fFilePath = os.path.join(self.fDirPath, self.fTypeOfCode_lower + "/combined_files")
        self.updateLatestFiles()
        self._readData()

    def _readData_PHITS(self, isCombined):
        def makeCombinedFiles():
            ResultDict = dict() # ResultDict[filename] =[(organID, equivDose, relativeError), ... (...)]
            # get data 
            for folderPath in self.fLatestFiles:
                # print("\n", folderPath)
                ResultDict[folderPath] = dict()
                # print(folderPath)
                # ----------------------------- neutrons ----------------------------- #
                # organID-WeigthSum.out...
                lst = natsort.natsorted([os.path.join(folderPath, _path) for _path in (os.listdir(folderPath)) if ((_path.endswith("WS.out")) or (\
                    _path.endswith("WS.out")))])
                # print(len(lst))
                for num in range(len(lst)): # 187 files
                    # print (f"{num}th iteration started |n")
                    f = open(lst[num], "r") # open file with that name
                    # print (lst[num], "is opend")
                    while True:
                        # read lines and split it by space so that I can get the list of each of words in lines
                        line = f.readline();

                        if not line: break; # if the line reaches to eof, then break
                        list_line = line.split()

                        if len(list_line)==0: # skip the \n lines
                            continue    
                        # neutron Parts
                        # fine line 
                        if (list_line[-1] == "r.err"): 
                            # num reg volume neutronDose r.err all r.err
                            list_line = f.readline().split()
                            cellID = list_line[1]
                            equivDose = float(list_line[3]); absDoseUncertainity = equivDose * float(list_line[4])
                            # check there a temp[cellID]
                            try: 
                                ResultDict[folderPath][cellID]
                            except:
                                ResultDict[folderPath][cellID] = [0., 0.] # temp[cellID] = [totalDose, relative Error]
                            ResultDict[folderPath][cellID][0] += equivDose;
                            ResultDict[folderPath][cellID][1] = math.sqrt(ResultDict[folderPath][cellID][1]**2 + absDoseUncertainity**2)
                            # print(ResultDict[folderPath][cellID][0], ResultDict[folderPath][cellID][1])
                    f.close()
                # ----------------------------- neutrons ----------------------------- #

                # ----------------------------- others ----------------------------- #
                lst = natsort.natsorted([os.path.join(folderPath, _path) \
                    for _path in (os.listdir(folderPath)) \
                        if (not (_path.endswith("WeightSum.out"))\
                        and not (_path.endswith("WS.out")))])[1:-1]
                for num in range(len(lst)): # 3 files
                    # print (f"{num}th iteration started |n")
                    filePath = lst[num]; GEM=False; HP=False; ALPHA=False
                    if "GEM" in filePath: GEM = True; HP=False; ALPHA=False
                    elif "HP" in filePath: GEM=False; HP=True; ALPHA=False
                    elif "ALPHA" in filePath: GEM=False; HP=False; ALPHA=True
                    f = open(filePath, "r") # open file with that name
                    # print (lst[num], "is opend")
                    while True:
                        # read lines and split it by space so that I can get the list of each of words in lines
                        line = f.readline();

                        if not line: break; # if the line reaches to eof, then break
                        list_line = line.split()

                        if len(list_line)==0: # skip the \n lines
                            continue    
                        # neutron Parts
                        # fine line 
                        if (list_line[-1] == "r.err"): 
                            # num reg volume neutronDose r.err all r.err
                            while True:
                                regNum = list_line[1];
                                if regNum.startswith("100"):
                                    break;
                                else:
                                    list_line = f.readline().split()
                            for i in range(187):
                                cellID = list_line[1];
                                # set equivDose and abs uncer. of it depending on the file type
                                if GEM:
                                    PP = float(list_line[3]); PPUncer = PP * float(list_line[4])
                                    PM = float(list_line[5]); PMUncer = PM * float(list_line[6])
                                    Elec = float(list_line[7]); ElecUncer = Elec * float(list_line[8])
                                    Photon = float(list_line[9]); PhotonUncer = Photon * float(list_line[10])
                                    equivDose = PP + PM + Elec + Photon
                                    absDoseUncertainity = math.sqrt(PPUncer**2 + PMUncer**2 + ElecUncer**2 + PhotonUncer**2)
                                elif HP:
                                    Proton = float(list_line[3]); ProtonUncer = Proton * float(list_line[4])
                                    PionP = float(list_line[5]); PionPUncer = PionP * float(list_line[6])
                                    PionM = float(list_line[7]); PionMUncer = PionM * float(list_line[8])
                                    equivDose = 2 * (Proton + PionP + PionM)
                                    absDoseUncertainity = math.sqrt(ProtonUncer**2 + PionPUncer**2 + PionMUncer**2)
                                elif ALPHA:
                                    Alpha = float(list_line[3]); AlphaUncer = Alpha * float(list_line[4])
                                    Nucleus = float(list_line[5]); NucleusUncer = Nucleus * float(list_line[6])
                                    equivDose = 20 * (Alpha + Nucleus)
                                    absDoseUncertainity = math.sqrt(AlphaUncer**2 + NucleusUncer**2)
                                # update result values
                                ResultDict[folderPath][cellID][0] += equivDose;
                                ResultDict[folderPath][cellID][1] = math.sqrt(ResultDict[folderPath][cellID][1]**2 + absDoseUncertainity**2)
                                list_line = f.readline().split()
                            # print(ResultDict[folderPath][cellID][0], ResultDict[folderPath][cellID][1])
                    f.close()
                # ----------------------------- others ----------------------------- #        
                
                ## set resultFolder
                today = date.today()
                ## Data file Dir
                dataFileDir = self.fDirPath
                ## output files
                phitsFolderDir = os.path.join(dataFileDir, "phits")
                ## resultFolderPath
                resultFolderPath = os.path.join(phitsFolderDir, f"combined_files/{today}")
                ## check dir 
                if not os.path.isdir(resultFolderPath):
                    os.makedirs(resultFolderPath)
                ## set the resultFileName
                dataType = folderPath.split("/")[-1]

                # for files in phitsDataFilesList:
                postureLst = ["op", "sp"]
                positionLst = [str(100 - pos*10) for pos in range(26)]
                positionLst = positionLst + ["7", "5", "3"]
                for posture in postureLst:
                    for position in positionLst:
                        text = f"{posture}_{position}"
                        if text == dataType:
                            resultFile = os.path.join(resultFolderPath, f"{posture}_{position}.out")

                # absError to relError
                for organID in ResultDict[folderPath].keys():
                    if ResultDict[folderPath][organID][0] != 0:
                        ResultDict[folderPath][organID][1] = ResultDict[folderPath][organID][1] / ResultDict[folderPath][organID][0]
                    else:
                        ResultDict[folderPath][organID][1] = 0
                        
                    
                g=open(resultFile, "w")
                g.write("organID \t equiv_Dose (Sv/source) \t relative error \n")
                for organID in ResultDict[folderPath].keys():
                    # '%-10s' % 'python'
                    sentence = '%-14s' % (str(organID)) + " | " + "{:.6e}".format(ResultDict[folderPath][organID][0])\
                        + "\t" +  "{:.6e}".format(ResultDict[folderPath][organID][1])
                    g.write(sentence + "\n")
                g.close()   

        if isCombined:
            makeCombinedFiles()
            
        self.fFilePath = os.path.join(self.fDirPath, self.fTypeOfCode_lower + "/combined_files")
        self.updateLatestFiles()
        self._readData()

    def calcEffDose(self, TWF: dict()):
        organIdToWf = TWF
        Result = dict()
        for file_name in self.fRawDataDic:
            # clear the totalEffdose list and organ Ed
            totalEffdose = [0., 0.]
            organEd = dict()
            # open the file
            with open(file_name) as f:
                # read line
                line = f.readline();
                
                while True:
                    # read the files until reach to the end of file
                    try:
                        line = f.readline().split();
                        # the line reaches to the first of organ data. 
                        # if line.startswith("100"):
                        if line[0] == "100":
                            # read all organs (187 organs)
                            for i in range(187):
                                # Get the organ ID, equivDose of it, and convert relative error to absolute error
                                organID = line[0];
                                # get equivalent Dose and calculate the abs uncertainity from the file data
                                equiv_dose = float(line[-2]); abs_error = float(line[-1]) * equiv_dose
                                # print(organID, equiv_dose, abs_error)
                                # find the TWF and multiply it to the equivDose, the organ is target organ

                                if organID in organIdToWf:
                                    # Set the orgnaEd[organId] to empty list
                                    organWeightingFactor = organIdToWf[organID][0]; organName = organIdToWf[organID][1]
                                    # organ Effective Dose dictionary, if if has not established, set the value to zero
                                    if organWeightingFactor == "r": ## organIdToWf[organID][0] = (WeightingFacotr, organName);   
                                        # if the organ is a remainder ,we should sum it up
                                        effDose = equiv_dose * (0.12 / 13.0); 
                                        effDoseAbsUncer = (abs_error * (0.12 / 13.0))
                                        # _individualEffDose = equiv_dose * 0.12; _individualEffDoseAbsUncer = (abs_error * 0.12)
                                        _individualEffDose = effDose; _individualEffDoseAbsUncer = effDoseAbsUncer
                                    else:
                                        # convert organWeightingFactor to float type from sting type
                                        organWeightingFactor = float(organWeightingFactor)
                                        # sum up the organ's effctive dose to total effective dose
                                        effDose = equiv_dose * organWeightingFactor
                                        effDoseAbsUncer = (abs_error * organWeightingFactor)
                                        _individualEffDose = effDose; _individualEffDoseAbsUncer = effDoseAbsUncer
                                    if np.isnan(effDose) or np.isnan(effDoseAbsUncer):
                                        pass
                                    else:
                                        totalEffdose[0] += effDose
                                        totalEffdose[1] = math.sqrt(totalEffdose[1]**2 + effDoseAbsUncer**2); 
                                    # calculate the effective doses of each organ
                                    # if not organEd[organName]:
                                    try:
                                        organEd[organName]
                                    except:
                                        organEd[organName] = [0. , 0.] # effDose and uncertainities for that organ
                                        organEd[organName][0] += _individualEffDose;
                                        organEd[organName][1] = math.sqrt(organEd[organName][1]**2 + _individualEffDoseAbsUncer**2); 
                                    else:
                                        organEd[organName][0] += _individualEffDose;
                                        organEd[organName][1] = math.sqrt(organEd[organName][1]**2 + _individualEffDoseAbsUncer**2); 
                                line = f.readline().split()
                    except:
                        break;
            calEffDose = totalEffdose[0]
            # print(calEffDose, type(calEffDose))
            if calEffDose != 0:
                calEffDoseRelError = totalEffdose[1]/totalEffdose[0]
            else:
                calEffDoseRelError = 0.0
            Result[file_name] = [calEffDose, calEffDoseRelError, organEd]

        self.fEffectiveDose = Result
        self.setData(self.fEffectiveDose)

    def writeResultEffDose(self):
        # write ResultEffDose to the file. 
        ## effDoseDic[file_name] = [totalEffdose[0], totalEffdose[1], organEd]
        numbering = 0;
        # resultfile = outfiledir + "result_effdose/" + f"{date}_effDose{numbering}.out"
        today = date.today()
        resultDir = os.path.join(self.fDirPath, f"resultEffDose/{today}")
        if not os.path.isdir(resultDir):
            os.makedirs(resultDir)
            
        resultFile = os.path.join(resultDir, f"{self.fTypeOfCode}-{numbering}.out")
        while os.path.isfile(resultFile):
            numbering +=1;
            resultFile = os.path.join(resultDir, f"{self.fTypeOfCode}-{numbering}.out")

        df = dict();

        g = open(resultFile, "w");
        df[self.fTypeOfCode] = dict()
        for key in self.fEffectiveDose.keys():
            name = str(key.split("/")[-1]).replace(".out","");
            df[self.fTypeOfCode][name+"Dose"] = self.fEffectiveDose[key][0]
            df[self.fTypeOfCode][name+"Error"] = self.fEffectiveDose[key][1]
            g.write("effective dose \t | relative error \n")
            sentence = '%-14s' % (str(key) +"\n \t" + "{:.6e}".format(self.fEffectiveDose[key][0]) + " Sv/source " + "{:.6e}".format(self.fEffectiveDose[key][1]) + "\n\n")
            g.write(sentence)
        # g.write("\n\n ---------------------------------------------------------------------------------------- \n")
        g.close()
        os.system(f"open {resultFile}")

    def writeOrganEffDose(self):

        # write ResultEffDose to the file. 
        ## effDoseDic[file_name] = [totalEffdose[0], totalEffdose[1], organEd]
        numbering = 0;
        # resultfile = outfiledir + "result_effdose/" + f"{date}_effDose{numbering}.out"
        today = date.today()

        resultDir = os.path.join(self.fDirPath, f"organEffDose/{today}")
        if not os.path.isdir(resultDir):
            os.makedirs(resultDir)


        mcCodeName = self.fTypeOfCode
        for key in self.fEffectiveDose.keys():

            folderPath = key
            dataType = folderPath.split("/")[-1]
            dataType = str(dataType).replace(".out", "")
            print(dataType)
            for posture in self.fPostureList:
                for position in self.fPositionList:
                    numbering = 0
                    string = f"{posture}_{position}"
                    if string == dataType:
                    # if string in folderPath:
                        resultFile = os.path.join(resultDir, f"{mcCodeName}_{posture}_{position}_{numbering}.out")
                        while os.path.isfile(resultFile):
                            numbering +=1;
                            resultFile = os.path.join(resultDir, f"{mcCodeName}_{posture}_{position}_{numbering}.out")

                        g = open(resultFile, "w");
                        sentenceFormat = "%-24s%-25s\n"
                        sentence = sentenceFormat % ("organName", "EffDose (Sv/source)")
                        g.write(sentence)
                        
                        sentenceFormat = "%-28s%-25s%-20s\n"
                        for organName in self.fEffectiveDose[key][2]:
                            sentence = sentenceFormat % (organName, \
                                "{:.6e}".format(self.fEffectiveDose[key][2][organName][0]), \
                                    "{:.6e}".format(self.fEffectiveDose[key][2][organName][1]))
                            g.write(sentence)
                        g.close()
                    
                    else:
                        pass

    def computeEDR(self):
        # variables
        activity = 3
        nps = 2.4 * 10**6
        mul = activity * nps
        mul = mul * 1.e6

        # fEffectiveDose -> {"fileName" : [effDose, effDoseError, organED]]}
        # fEDRs -> {"fileName": [EDRs, EDRsError]}
        # organContr will compute each organ's contribution of EDRs

        for items in self.fEffectiveDose.items():
            fileName = items[0]
            data = items[1]

            case = fileName.split("/")[-1].replace(".out", "");

            effDose = data[0]
            effDoseError = data[1]

            self.fEDRs[case] = [effDose * mul, effDoseError]
            self.setDataFrame_EDRs()
            self.sortEDRs()
            # unit = Sv/second

    def convertOrganEffDoseContr(self):
        activity = 3
        nps = 2.4 * 10**6
        mul = activity * nps
        mul = mul * 1.e6

        for items in self.fEffectiveDose.items():
            fileName = items[0]
            data = items[1]
            
            case = fileName.split("/")[-1].replace(".out", "");
            self.fOrganContr[case] = dict()

            # organED = dict {"organName" : [contribution, AbsError**2]}
            organED = data[2]
            # print(organED.items())
            for organName, valueLst in organED.items():
                organEDR = valueLst[0] * mul
                organEDRError = valueLst[1] * mul
                # unit = Sv/second
                self.fOrganContr[case][organName] = [organEDR, organEDRError]
        self.setDataFram_OrganContr()
        self.sortORganEDRsContr()

    def sortEDRs(self):
        ### self.fEDRs["op_100"] = [effDose, effDoseError]
        indexlst = ["op_100", "op_90", "op_80", "op_70", "op_60", "op_50", "op_40", "op_30", "op_20", "op_10", "op_7", "op_5", "op_3", "op_0", \
                "op_-10", "op_-20", "op_-30", "op_-40", "op_-50", "op_-60", "op_-70", "op_-80", "op_-90", "op_-100", "op_-110", "op_-120", "op_-130", "op_-140", "op_-150", \
                "sp_100", "sp_90", "sp_80", "sp_70", "sp_60", "sp_50", "sp_40", "sp_30", "sp_20", "sp_10", "sp_7", "sp_5", "sp_3", "sp_0", \
                "sp_-10", "sp_-20", "sp_-30", "sp_-40", "sp_-50", "sp_-60", "sp_-70", "sp_-80", "sp_-90", "sp_-100", "sp_-110", "sp_-120", "sp_-130", "sp_-140", "sp_-150", \
                "sq_100", "sq_90", "sq_80", "sq_70", "sq_60", "sq_50", "sq_40", "sq_30", "sq_20", "sq_10", "sq_7", "sq_5", "sq_3", "sq_0", \
                "sq_-10", "sq_-20", "sq_-30", "sq_-40", "sq_-50", "sq_-60", "sq_-70", "sq_-80", "sq_-90", "sq_-100", "sq_-110", "sq_-120", "sq_-130", "sq_-140", "sq_-150", \
                "bd_100", "bd_90", "bd_80", "bd_70", "bd_60", "bd_50", "bd_40", "bd_30", "bd_20", "bd_10", "bd_7", "bd_5", "bd_3", "bd_0", \
                "bd_-10", "bd_-20", "bd_-30", "bd_-40", "bd_-50", "bd_-60", "bd_-70", "bd_-80", "bd_-90", "bd_-100", "bd_-110", "bd_-120", "bd_-130", "bd_-140", "bd_-150"]

        self.fDF_EDRs = self.fDF_EDRs.reindex(index=indexlst, copy=True)

        self.fDF_EDRs.insert(2, "type", ["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
                          "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
                          "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
                          "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
        df = self.fDF_EDRs.groupby(self.fDF_EDRs.type)

        df_op = df.get_group("op")
        df_sp = df.get_group("sp")
        df_sq = df.get_group("sq")
        df_bd = df.get_group("bd")

        self.fDF_EDRs = dict()
        self.fDF_EDRs["OP"] = df_op
        self.fDF_EDRs["SP"] = df_sp
        self.fDF_EDRs["SQ"] = df_sq
        self.fDF_EDRs["BD"] = df_bd

        for key in self.fDF_EDRs.keys():
            self.fDF_EDRs[key].drop(columns="type", inplace=True)

    def sortORganEDRsContr(self):
        ### self.fEDRs["op_100"] = [effDose, effDoseError]
        indexlst = ["op_100", "op_90", "op_80", "op_70", "op_60", "op_50", "op_40", "op_30", "op_20", "op_10", "op_7", "op_5", "op_3", "op_0", \
                "op_-10", "op_-20", "op_-30", "op_-40", "op_-50", "op_-60", "op_-70", "op_-80", "op_-90", "op_-100", "op_-110", "op_-120", "op_-130", "op_-140", "op_-150", \
                "sp_100", "sp_90", "sp_80", "sp_70", "sp_60", "sp_50", "sp_40", "sp_30", "sp_20", "sp_10", "sp_7", "sp_5", "sp_3", "sp_0", \
                "sp_-10", "sp_-20", "sp_-30", "sp_-40", "sp_-50", "sp_-60", "sp_-70", "sp_-80", "sp_-90", "sp_-100", "sp_-110", "sp_-120", "sp_-130", "sp_-140", "sp_-150", \
                "sq_100", "sq_90", "sq_80", "sq_70", "sq_60", "sq_50", "sq_40", "sq_30", "sq_20", "sq_10", "sq_7", "sq_5", "sq_3", "sq_0", \
                "sq_-10", "sq_-20", "sq_-30", "sq_-40", "sq_-50", "sq_-60", "sq_-70", "sq_-80", "sq_-90", "sq_-100", "sq_-110", "sq_-120", "sq_-130", "sq_-140", "sq_-150", \
                "bd_100", "bd_90", "bd_80", "bd_70", "bd_60", "bd_50", "bd_40", "bd_30", "bd_20", "bd_10", "bd_7", "bd_5", "bd_3", "bd_0", \
                "bd_-10", "bd_-20", "bd_-30", "bd_-40", "bd_-50", "bd_-60", "bd_-70", "bd_-80", "bd_-90", "bd_-100", "bd_-110", "bd_-120", "bd_-130", "bd_-140", "bd_-150"]

        self.fDF_OrganContr = self.fDF_OrganContr.reindex(index=indexlst, copy=True)

        self.fDF_OrganContr.insert(0, "type", ["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
                          "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
                          "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
                          "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
        df = self.fDF_OrganContr.groupby(self.fDF_OrganContr.type)

        df_op = df.get_group("op")
        df_sp = df.get_group("sp")
        df_sq = df.get_group("sq")
        df_bd = df.get_group("bd")

        self.fDF_OrganContr = dict()
        self.fDF_OrganContr["OP"] = df_op
        self.fDF_OrganContr["SP"] = df_sp
        self.fDF_OrganContr["SQ"] = df_sq
        self.fDF_OrganContr["BD"] = df_bd

        for key in self.fDF_OrganContr.keys():
            self.fDF_OrganContr[key].drop(columns="type", inplace=True)

    # --- Plotting region ---
    def init_plot(self, imageFolderPath: str, figsize = (12,8)):
        figure = plt.figure(figsize=figsize)
        self.figure = figure
        self.fmt = ["o", "^", "s"]

        ax = figure.add_subplot(111)
        self.ax = ax

        self.fOptions = {"legendFontSize": 20,
                         "xLabeltxt": "Position (cm)",
                         "yLabeltxt": "Effective Dose Rate (\u03BCSv/s)",
                         "labelFontSize": 15,
                         "xticksFontSize": 15,
                         "xScale": "linear",
                         "yScale": "linear",
                         "title": "Effective Dose Rate for Operating Worker",
                         "titleFontSize": 20,
                         "grid": {"alpha": 0.7, "color": "gray", "linestyle": "--"},
                         "xLim": [-160, 110],
                         "yLim": [5.e-07, 1]
                         }

        self.today = str(date.today())
        self.imageFolerDir = os.path.join(imageFolderPath, self.today)
        self.saveForlderPath = ""
    
    def refreshFigure(self):
        """ This functions refresh the setups and redraw"""
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    # --- Set methods ---
    def setX(self, X: list):
        """ inputs lists for X """
        self.x = X

    def setY(self, Y: list):
        """ inputs lists for Y """
        self.y = Y

    def setYError(self, YError: list):
        """ inputs lists for Absoulte Error of Y """
        self.yError = YError

    def setOptions(self, OptionDict: dict):
        """ Set Options as input Option Dictionary """
        self.fOptions = OptionDict

    def replaceOption(self, optionNameLst: list, valueLst: list):
        """ Replace some of options as the input values """
        for optionName, value in zip(optionNameLst, valueLst):
            self.fOptions[optionName] = value

    def applyOptions(self):
        """ Apply the options for plotting with self.fOptions """
        # plt.legend("GEANT4")
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
        self.refreshFigure()

    def setToLinX(self):
        """ Set X Scale to Linear"""
        self.ax.set_xscale("linear")
        self.refreshFigure()

    def setToLogX(self):
        """ Set X Scale to Log """
        self.ax.set_xscale("log")
        self.refreshFigure()

    def setToLinY(self):
        """ Set Y Scale to Linear """
        self.ax.set_yscale("linear")
        self.refreshFigure()

    def setToLogY(self):
        """ Set Y Scale to Log """
        self.ax.set_yscale("log")
        self.refreshFigure()

    def getX(self):
        return self.fPositionList
    
    def getY(self, posture = "OP"):
        return self.getDF_EDRs(posture)["EDR"].values
        
    def getYErr(self, posture = "OP"):
        return self.getDF_EDRs(posture)["EDR_RelErr"].values * self.getDF_EDRs(posture)["EDR"].values

    def getOptions(self):
        return self.fOptions

    # --- Draw operating, supporting, bending, and squatting --- 
    def _setSaveFolder(self, folderPath):
        """ sub function for setting save folder path and make directories """
        newFolderPath = os.path.join(self.imageFolerDir, folderPath)
        self.saveForlderPath = newFolderPath

        try: 
            os.makedirs(self.saveForlderPath)
        except:
            pass

    def saveFig(self, linlog = True, posture="OP", position="GROUND"):
        """ save picture under the saveFolder Path at: FolderPath/Trendlines/<Posture> """
        posture = posture.upper()
        codeName = self.fTypeOfCode
        folderPath = os.path.join("Trendlines", posture)
        self._setSaveFolder(folderPath=folderPath)

        if linlog:
            self.setToLinY()
            self.figure.savefig( os.path.join(self.saveForlderPath, f"{codeName}_{posture}_{position}_LIN.svg") )
            self.setToLogY()
            self.figure.savefig( os.path.join(self.saveForlderPath, f"{codeName}_{posture}_{position}_LOG.svg") )
        else:
            self.figure.savefig( os.path.join(self.saveForlderPath, f"{codeName}_{posture}_{position}_LIN.svg") )
        plt.close()

    def drawErrorBar(self, posture, fmt = "o", field = "Total", xscale = "linear", yscale = "linear", isSave = False):
        if field == "Total":
            self.refreshFigure()

        elif field == "Ground":
            pass

        elif field == "Underground":
            pass

        codeName = self.fTypeOfCode
        # x = self.getX()
        # y = self.getY(posture = posture)
        # yerr = self.getYErr(posture=posture)
        self.setX( self.getX() )
        self.setY( self.getY(posture = posture) )
        self.setYError( self.getYErr(posture=posture) )
        # print(codeName, posture, np.isnan(self.y)[0])
        if np.isnan(self.y)[0] == True:
            print("No data!")
            return 0
        
        self.ax.errorbar(x = self.x, y = self.y, yerr = self.yError, label = codeName, fmt = fmt)
        self.replaceOption(["title"], [f"Effective Dose Rate for {posture} Worker"])
        self.applyOptions()

        if xscale == "linear":
            self.setToLinX()
        else:
            self.setToLogX()
        
        if yscale == "linear":
            self.setToLinY()
        else:
            self.setToLogY()
        self.refreshFigure()

        if isSave:
            if xscale == "linear" and yscale == "log":
                self.saveFig(linlog=True, posture = posture, position = field)
            else:
                self.saveFig(linlog=False, posture = posture, position = field)
        else:
            pass
        
    def drawBasicFigures(self):
        for posture in self.getPostureList():
            self.drawErrorBar(posture=posture, yscale="log", isSave=True)
            self.ax.cla()

    def exportPlotting(self, posture):
        x = self.getX()
        y = self.getY(posture=posture)
        yerr = self.getYErr(posture=posture)
        codeName = self.getCodeName()
        return (x, y, yerr, codeName)

    
    ### fitting region
    def getCoef(self, posture, data):
        self.fCoef[posture] = data



class SupplementaryData(Data):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # --- * --- Data Variables
        self.initialize()

    def initialize(self):
        print(self.fFilePath)
        with open(self.fFilePath, "r") as ef:
            while True:
                try:
                    # read each line of file
                    line = ef.readline().strip().split(",")
                    # set variables
                    organName = line[0]; WeightingFactor = line[1]; organID_list = line[2:];
                    for organID in organID_list:
                        self.fData[organID] = (WeightingFactor, organName);
                except:
                    break;
    
    def getData(self):
        return self.fData      
