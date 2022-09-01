import math
import os
import natsort
from collections import OrderedDict
from datetime import date
import numpy as np
import pandas as pd
import pickle

def getMcnp6FileFolders():
    ## Data file Dir
    dataFileDir = "../../data_source/data/sonde/"
    ## output files
    mcnp6FolderDir = os.path.join(dataFileDir, "mcnp6")
    ## Get the path fo latest files' directory
    mcnp6RecentFolderPath =  os.path.join(mcnp6FolderDir ,\
        natsort.natsorted([folder for folder in (os.listdir(mcnp6FolderDir)) if folder.startswith("202")])[-1])
    ## get the latest folders, "../ ... /equip", "../ ... /s_1m"
    mcnp6RecentFolderPath = [os.path.join(mcnp6RecentFolderPath, _path) \
        for _path in natsort.natsorted(os.listdir(mcnp6RecentFolderPath), reverse=True)]
    
    return mcnp6RecentFolderPath

def makeCombindedFilesMcnp6(mcnp6DataFilesList: list):
    nps1e6 = "1000000"; nps1e7 = "10000000"; nps2e5 = "200000"; nps3e7 = "30000000"; nps5e5 = '500000'; nps25e5 = "250000"; nps5e6 ="5000000"; nps1e5 = "100000"
    npsList = [nps1e6, nps1e7, nps2e5, nps3e7, nps5e5, nps25e5, nps5e6, nps1e5]
    ResultDict = OrderedDict() # ResultDict[filename] =[(organID, equivDose, relativeError), ... (...)]
    # get data 
    
    for folderPath in mcnp6DataFilesList:
        
        ResultDict[folderPath] = OrderedDict()
        # equip..0 - equip..9
        lst = [os.path.join(folderPath, _path) for _path in (os.listdir(folderPath))]
        # print(folderPath)
        # print(lst)
        for num in range(len(lst)):
            temp = OrderedDict() # temp[cellID] = [totalDose, relative Error]
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

        
        # for key in ResultDict[folderPath].keys():
        #     print(ResultDict[folderPath][key])
        
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
        dataFileDir = "../../data_source/data/sonde/"
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

def getPhitsFileFolders():
    ## Data file Dir
    dataFileDir = "../../data_source/data/sonde/"
    ## output files
    phitsFolderDir = os.path.join(dataFileDir, "phits")
    ## Get the path fo latest files' directory
    phitsRecentFolderPath =  os.path.join(phitsFolderDir ,\
        natsort.natsorted([folder for folder in (os.listdir(phitsFolderDir)) if folder.startswith("202")])[-1])
    ## get the latest folders, "../ ... /equip", "../ ... /s_1m"
    phitsRecentFolderPath = [os.path.join(phitsRecentFolderPath, _path) \
        for _path in natsort.natsorted(os.listdir(phitsRecentFolderPath), reverse=True)]
    
    return phitsRecentFolderPath
    
def makeCombinedFilesPhits(phitsDataFilesList: list):
    ResultDict = OrderedDict() # ResultDict[filename] =[(organID, equivDose, relativeError), ... (...)]
    # get data 
    for folderPath in phitsDataFilesList:
        # print("\n", folderPath)
        ResultDict[folderPath] = OrderedDict()
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
        dataFileDir = "../../data_source/data/sonde/"
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
    
def getFiles():
    # Set Variables
    ## data File Directory
    dataFileDir = "../../data_source/data/sonde/"

    ## output files or handled files directory 
    geant4FolderDir = os.path.join(dataFileDir, "geant4"); 
    mcnp6FolderDir = os.path.join(dataFileDir, "mcnp6/combined_files")
    phitsFolderDir = os.path.join(dataFileDir, "phits/combined_files")

    ## get the latest files' directory
    geant4RecentFolderPath =  os.path.join(geant4FolderDir ,natsort.natsorted([folder for folder in (os.listdir(geant4FolderDir)) if folder.startswith("202")])[-1])
    mcnp6RecentFolderPath =  os.path.join(mcnp6FolderDir ,natsort.natsorted([folder for folder in (os.listdir(mcnp6FolderDir)) if folder.startswith("202")])[-1])
    phitsRecentFolderPath =  os.path.join(phitsFolderDir ,natsort.natsorted([folder for folder in (os.listdir(phitsFolderDir)) if folder.startswith("202")])[-1])

    ## print out
    # print(f" Geant4 data path : {geant4RecentFolderPath}")
    # print(f" MCNP6  data path : {mcnp6RecentFolderPath}")

    ## get the latest files
    geant4LatestFiles = [os.path.join(geant4RecentFolderPath, _path) \
        for _path in natsort.natsorted(os.listdir(geant4RecentFolderPath), reverse=True)]
    mcnp6LatestFiles = [os.path.join(mcnp6RecentFolderPath, _path) \
        for _path in natsort.natsorted(os.listdir(mcnp6RecentFolderPath), reverse=True)]
    phitsLatestFiles = [os.path.join(phitsRecentFolderPath, _path) \
        for _path in natsort.natsorted(os.listdir(phitsRecentFolderPath), reverse=True)]

    ## print out 
    # print (geant4LatestFiles)
    # print (mcnp6LatestFiles)
    # print (phitsLatestFiles)

    return geant4LatestFiles, mcnp6LatestFiles, phitsLatestFiles

def getTWF():
    ### Get the tissue weighting factors of each organ ###
    organIdToWf = dict()
    # effetive_tissue_IDs data_source/load_data
    eff_path ="../../data_source/load_data/eff_tissue_IDs.inp"
    # Read files, foramt = organ name, WF, orID1, orID2, ... , orIDn
    with open(eff_path, "r") as ef:
        while True:
            try:
                # read each line of file
                line = ef.readline().strip().split(",")
                # set variables
                organName = line[0]; WeightingFactor = line[1]; organID_list = line[2:];
                for organID in organID_list:
                    organIdToWf[organID] = (WeightingFactor, organName);
            except:
                break;
    return organIdToWf;
                     
def calEffectiveDose(dataFile, organIdToWf: dict()) -> dict():
    # dataFile => list of file names of data, getFiles ordered
    # / orgnaIdToWf => dictionary of getTWF, organIdToWf[organID] = (WeightingFactor, organName);   
    # set Variables
    Result = OrderedDict()
    for file_name in dataFile:
        # clear the totalEffdose list and organ Ed
        totalEffdose = [0., 0.]
        organEd = OrderedDict()
        # open the file
        # print(file_name)
        # print(dataFile)
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
                                # else:
                                    # print("abs")
                                # else:
                                #     organEd[organName][0] += _individualEffDose;
                                #     organEd[organName][1] = math.sqrt(organEd[organName][1]**2 + _individualEffDoseAbsUncer**2); 
                                    # print(organEd[organName])
                            line = f.readline().split()
                except:
                    break;
        # result[filename] = [Effdose, Relative Uncertainity, Effective dose dictionary of organs]
        # print(totalEffdose[0], totalEffdose[1], organEd)
        calEffDose = totalEffdose[0]
        if not calEffDose == 0 or calEffDose == "-nan":
            calEffDoseRelError = totalEffdose[1]/totalEffdose[0]
        else:
            calEffDoseRelError = 0.0
        # Result[file_name] = [totalEffdose[0], totalEffdose[1]/totalEffdose[0], organEd]
        Result[file_name] = [calEffDose, calEffDoseRelError, organEd]
        # print(organEd.items())
    return Result

def writeResultEffDose(*args):
    # write ResultEffDose to the file. 
    ## effDoseDic[file_name] = [totalEffdose[0], totalEffdose[1], organEd]
    numbering = 0;
    # resultfile = outfiledir + "result_effdose/" + f"{date}_effDose{numbering}.out"
    today = date.today()
    resultDir = f"../../data_source/data/sonde/resultEffDose/{today}"
    if not os.path.isdir(resultDir):
        os.makedirs(resultDir)
        # print("error"
        
    resultFile = os.path.join(resultDir, f"{today}-{numbering}.out")
    while os.path.isfile(resultFile):
        numbering +=1;
        resultFile = os.path.join(resultDir, f"{today}-{numbering}.out")
    df = OrderedDict();

    g = open(resultFile, "w");
    strings= ["GEANT4", "MCNP6", "PHTIS"]
    indexing = 0
    for arg in args:
        df[strings[indexing]] = OrderedDict()
        for key in arg.keys():
            name = str(key.split("/")[-1]).replace(".out","");
            df[strings[indexing]][name+"Dose"] = arg[key][0]
            df[strings[indexing]][name+"Error"] = arg[key][1]
            g.write("effective dose \t | relative error \n")
            sentence = '%-14s' % (str(key) +"\n \t" + "{:.6e}".format(arg[key][0]) + " Sv/source " + "{:.6e}".format(arg[key][1]) + "\n\n")
            g.write(sentence)
        g.write("\n\n ---------------------------------------------------------------------------------------- \n")
    g.close()
    df = pd.DataFrame(df)
    # print(df)
    os.system(f"open {resultFile}")
    
def writeOrganEffectiveDose(*args):
    
    # write ResultEffDose to the file. 
    ## effDoseDic[file_name] = [totalEffdose[0], totalEffdose[1], organEd]
    numbering = 0;
    # resultfile = outfiledir + "result_effdose/" + f"{date}_effDose{numbering}.out"
    today = date.today()
    resultDir = f"../../data_source/data/sonde/organEffDose/{today}"
    if not os.path.isdir(resultDir):
        os.makedirs(resultDir)
        # print("error"
    
    strings= ["GEANT4", "MCNP6", "PHITS"]; strIndex = 0
    postureList = ["op", "sp", "sq", "bd"]
    # positionList = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0, -30, -40, -50, -60, -70, -80, -90, -100, -110, -120, -130, -140, -150]
    positionList = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 7, 5, 3, 0, -30, -40, -50, -60, -70, -80, -90, -100, -110, -120, -130, -140, -150]
    for arg in args:
        mcCodeName = strings[strIndex]; strIndex+=1
        for key in arg.keys():

            folderPath = key
            dataType = folderPath.split("/")[-1]
            dataType = str(dataType).replace(".out", "")
            
            print(mcCodeName)
            print(dataType)


            for posture in postureList:
                for position in positionList:
                    numbering = 0
                    string = str(posture)+"_"+str(position)
                    string = f"{posture}_{position}"
                    if string == dataType:
                    # if string in folderPath:
                        resultFile = os.path.join(resultDir, f"{today}_{mcCodeName}_{posture}_{position}_{numbering}.out")
                        while os.path.isfile(resultFile):
                            numbering +=1;
                            resultFile = os.path.join(resultDir, f"{today}_{mcCodeName}_{posture}_{position}_{numbering}.out")
                        g = open(resultFile, "w");
                        sentenceFormat = "%-24s%-25s\n"
                        sentence = sentenceFormat % ("organName", "EffDose (Sv/source)")
                        g.write(sentence)
                        
                        sentenceFormat = "%-28s%-25s%-20s\n"
                        for organName in arg[key][2]:
                            sentence = sentenceFormat % (organName, \
                                "{:.6e}".format(arg[key][2][organName][0]), "{:.6e}".format(arg[key][2][organName][1]))
                            g.write(sentence)
                        g.close()
                    
                    else:
                        # resultFile = os.path.join(resultDir, f"{today}_{numbering}.out")
                        pass                    
                    
            
def getAnnualDose(*args):
    numbering = 0;
    today = date.today()
    resultDir = f"../../data_source/data/sonde/AnnualDose/{today}"
    if not os.path.isdir(resultDir):
        os.makedirs(resultDir)
    
    resultFile = os.path.join(resultDir, f"{today}_{numbering}.out")
    while os.path.isfile(resultFile):
        numbering +=1;
        resultFile = os.path.join(resultDir, f"{today}_{numbering}.out")
    
    data = OrderedDict(); dataRel = OrderedDict();

    data["GEANT4"] = dict(); 
    data["MCNP6"] = dict(); 
    data["PHITS"] = dict()

    dataRel["GEANT4"] = dict()
    dataRel["MCNP6"] = dict()
    dataRel["PHITS"] = dict()
    columns = []

    ### variables 
    activity = 3
    nps = 2.4 * 10**6
    mul = activity * nps

    ### calculations
    for arg in args:
        ### select files geant4Result, mcnp6Result, phitsResult
        for key in arg.keys():
            keylst = key.split("/");
            codes = keylst[5]; case = keylst[-1].replace(".out", "");

            if (not "CellFlux" in case):
                if codes == "geant4":
                    data["GEANT4"][case] = arg[key][0] * mul
                    dataRel["GEANT4"][case]= arg[key][1]

                elif codes == "mcnp6":
                    data["MCNP6"][case] = arg[key][0] * mul
                    dataRel["MCNP6"][case]= arg[key][1]

                elif codes == "phits":
                    data["PHITS"][case] = arg[key][0] * mul
                    dataRel["PHITS"][case]= arg[key][1]
                    
                columns.append(case)
    
    print(data)
    indexlst = ["op_100", "op_90", "op_80", "op_70", "op_60", "op_50", "op_40", "op_30", "op_20", "op_10", "op_7", "op_5", "op_3", "op_0", \
                "op_-10", "op_-20", "op_-30", "op_-40", "op_-50", "op_-60", "op_-70", "op_-80", "op_-90", "op_-100", "op_-110", "op_-120", "op_-130", "op_-140", "op_-150", \
                "sp_100", "sp_90", "sp_80", "sp_70", "sp_60", "sp_50", "sp_40", "sp_30", "sp_20", "sp_10", "sp_7", "sp_5", "sp_3", "sp_0", \
                "sp_-10", "sp_-20", "sp_-30", "sp_-40", "sp_-50", "sp_-60", "sp_-70", "sp_-80", "sp_-90", "sp_-100", "sp_-110", "sp_-120", "sp_-130", "sp_-140", "sp_-150", \
                "sq_100", "sq_90", "sq_80", "sq_70", "sq_60", "sq_50", "sq_40", "sq_30", "sq_20", "sq_10", "sq_7", "sq_5", "sq_3", "sq_0", \
                "sq_-10", "sq_-20", "sq_-30", "sq_-40", "sq_-50", "sq_-60", "sq_-70", "sq_-80", "sq_-90", "sq_-100", "sq_-110", "sq_-120", "sq_-130", "sq_-140", "sq_-150", \
                "bd_100", "bd_90", "bd_80", "bd_70", "bd_60", "bd_50", "bd_40", "bd_30", "bd_20", "bd_10", "bd_7", "bd_5", "bd_3", "bd_0", \
                "bd_-10", "bd_-20", "bd_-30", "bd_-40", "bd_-50", "bd_-60", "bd_-70", "bd_-80", "bd_-90", "bd_-100", "bd_-110", "bd_-120", "bd_-130", "bd_-140", "bd_-150"]
    # indexlst = ["op_100", "op_90", "op_80", "op_70", "op_60", "op_50", "op_40", "op_30", "op_20", "op_10", "op_0", \
    #             "op_-10", "op_-20", "op_-30", "op_-40", "op_-50", "op_-60", "op_-70", "op_-80", "op_-90", "op_-100", "op_-110", "op_-120", "op_-130", "op_-140", "op_-150", \
    #             "sp_100", "sp_90", "sp_80", "sp_70", "sp_60", "sp_50", "sp_40", "sp_30", "sp_20", "sp_10", "sp_0", \
    #             "sp_-10", "sp_-20", "sp_-30", "sp_-40", "sp_-50", "sp_-60", "sp_-70", "sp_-80", "sp_-90", "sp_-100", "sp_-110", "sp_-120", "sp_-130", "sp_-140", "sp_-150", \
    #             "sq_100", "sq_90", "sq_80", "sq_70", "sq_60", "sq_50", "sq_40", "sq_30", "sq_20", "sq_10", "sq_0", \
    #             "sq_-10", "sq_-20", "sq_-30", "sq_-40", "sq_-50", "sq_-60", "sq_-70", "sq_-80", "sq_-90", "sq_-100", "sq_-110", "sq_-120", "sq_-130", "sq_-140", "sq_-150", \
    #             "bd_100", "bd_90", "bd_80", "bd_70", "bd_60", "bd_50", "bd_40", "bd_30", "bd_20", "bd_10", "bd_0", \
    #             "bd_-10", "bd_-20", "bd_-30", "bd_-40", "bd_-50", "bd_-60", "bd_-70", "bd_-80", "bd_-90", "bd_-100", "bd_-110", "bd_-120", "bd_-130", "bd_-140", "bd_-150"]


    df = pd.DataFrame(data)
    df = pd.DataFrame(data, index=columns[:116])
    df = df.reindex(index=indexlst)
    # df.insert(1, "type",   ["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
    #                         "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
    #                         "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
    #                         "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
    df.insert(1, "type", ["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
                          "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
                          "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
                          "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
    df_group = df.groupby(df.type)
    
    dfRel = pd.DataFrame(dataRel, index=columns[:116])
    dfRel = dfRel.reindex(index=indexlst)
    # dfRel.insert(1, "type",["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
    #                         "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
    #                         "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
    #                         "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
    dfRel.insert(1, "type", ["op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", "op", \
                          "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", "sp", \
                          "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", "sq", \
                          "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd", "bd"])
    
    dfRel_group = dfRel.groupby(dfRel.type)

    
    df1 = df_group.get_group("op")
    df1.rename(columns={
        "GEANT4" : "GEANT4_op",
        "MCNP6" : "MCNP6_op",
        "PHITS" : "PHITS_op"
                        }, inplace=True)
    df1.drop(["type"], axis=1, inplace=True)
    df1.index = indexlst[:29]
    
    df2 = df_group.get_group("sp")
    df2.rename(columns={
        "GEANT4" : "GEANT4_sp",
        "MCNP6" : "MCNP6_sp",
        "PHITS" : "PHITS_sp"
                        }, inplace=True)
    df2.index = indexlst[:29]
    df2.drop(["type"], axis=1, inplace=True)

    df3 = df_group.get_group("sq")
    df3.rename(columns={
        "GEANT4" : "GEANT4_sq",
        "MCNP6" : "MCNP6_sq",
        "PHITS" : "PHITS_sq"
                        }, inplace=True)
    df3.index = indexlst[:29]
    df3.drop(["type"], axis=1, inplace=True)

    df4 = df_group.get_group("bd")
    df4.rename(columns={
        "GEANT4" : "GEANT4_bd",
        "MCNP6" : "MCNP6_bd",
        "PHITS" : "PHITS_bd"
                        }, inplace=True)
    df4.index = indexlst[:29]
    df4.drop(["type"], axis=1, inplace=True)
    
    dfRel1 = dfRel_group.get_group("op")
    dfRel1.rename(columns={
        "GEANT4" : "GEANT4_ERROR_op",
        "MCNP6" : "MCNP6_ERROR_op",
        "PHITS" : "PHITS_ERROR_op"
                            }, inplace=True)
    dfRel1.drop(["type"], axis=1, inplace=True)
    dfRel1.index = indexlst[:29]
    
    dfRel2 = dfRel_group.get_group("sp")
    dfRel2.rename(columns={
        "GEANT4" : "GEANT4_ERROR_sp",
        "MCNP6" : "MCNP6_ERROR_sp",
        "PHITS" : "PHITS_ERROR_sp"
                            }, inplace=True)
    dfRel2.index = indexlst[:29]
    dfRel2.drop(["type"], axis=1, inplace=True)

    dfRel3 = dfRel_group.get_group("sq")
    dfRel3.rename(columns={
        "GEANT4" : "GEANT4_ERROR_sq",
        "MCNP6" : "MCNP6_ERROR_sq",
        "PHITS" : "PHITS_ERROR_sq"
                            }, inplace=True)
    dfRel3.index = indexlst[:29]
    dfRel3.drop(["type"], axis=1, inplace=True)

    dfRel4 = dfRel_group.get_group("bd")
    dfRel4.rename(columns={
        "GEANT4" : "GEANT4_ERROR_bd",
        "MCNP6" : "MCNP6_ERROR_bd",
        "PHITS" : "PHITS_ERROR_bd"
                            }, inplace=True)
    dfRel4.index = indexlst[:29]
    dfRel4.drop(["type"], axis=1, inplace=True)
    
    df = pd.concat([df1, df2, df3, df4], axis=1)
    # print(df)
    # df.index = ["100cm", "90cm", "80cm", "70cm", "60cm", "50cm", "40cm", "30cm", "20cm", "10cm", "0cm", \
    df.index = ["100cm", "90cm", "80cm", "70cm", "60cm", "50cm", "40cm", "30cm", "20cm", "10cm", "7cm", "5cm", "3cm", "0cm", \
                "-10cm", "-20cm", "-30cm", "-40cm", "-50cm", "-60cm", "-70cm", "-80cm", "-90cm", "-100cm", "-110cm", "-120cm", "-130cm", "-140cm", "-150cm"]
    # print(df)
    
    dfRel = pd.concat([dfRel1, dfRel2, dfRel3, dfRel4], axis=1)
    # dfRel.index = ["100cm", "90cm", "80cm", "70cm", "60cm", "50cm", "40cm", "30cm", "20cm", "10cm", "0cm", \
    dfRel.index = ["100cm", "90cm", "80cm", "70cm", "60cm", "50cm", "40cm", "30cm", "20cm", "10cm", "7cm", "5cm", "3cm", "0cm", \
                   "-10cm", "-20cm", "-30cm", "-40cm", "-50cm", "-60cm", "-70cm", "-80cm", "-90cm", "-100cm", "-110cm", "-120cm", "-130cm", "-140cm", "-150cm"]
    # print(dfRel)
    
    dfResult = pd.concat([df, dfRel], axis=1)
    # dfResult.insert(0, "position", [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0, \
    dfResult.insert(0, "position", [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 7, 5, 3, 0, \
                                    -10, -20, -30, -40, -50, -60, -70, -80, -90, -100, -110, -120, -130, -140, -150])

    # resultlstOp = ["position", "GEANT4_op", "GEANT4_ERROR_op", "MCNP6_op", "MCNP6_ERROR_op", "PHITS_op", "PHITS_ERROR_op"]
    # resultlstSp = ["position", "GEANT4_sp", "GEANT4_ERROR_sp", "MCNP6_sp", "MCNP6_ERROR_sp", "PHITS_sp", "PHITS_ERROR_sp"]
    # dfResultOp = dfResult[resultlstOp]
    # dfResultSp = dfResult[resultlstSp]

    # print(dfResult)
    # print(dfResultOp)
    # print(dfResultSp)
    
    dfResult.to_csv(resultFile)
    with open(f"../../pickle/AnnualDose_{today}", "wb") as file:
        pickle.dump(dfResult, file)