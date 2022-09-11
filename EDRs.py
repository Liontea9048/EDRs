from EDRs_pkg._subEquivToEffective import *

if __name__ == "__main__":
    # handle the mcnp6 files into combined files
    # mcnp6FolderPath = getMcnp6FileFolders()
    # makeCombindedFilesMcnp6(mcnp6FolderPath)
    
    # handle the phits files into combined files
    # phitsFolderPath =getPhitsFileFolders()
    # makeCombinedFilesPhits(phitsFolderPath)
    
    # # calculate effective dose
    geant4Data, mcnp6Data, phitsData = getFiles()
    organTWF = getTWF()
    geant4Result = calEffectiveDose(geant4Data, organTWF)
    # mcnp6Result = calEffectiveDose(mcnp6Data, organTWF)
    # phitsResult = calEffectiveDose(phitsData, organTWF)

    # writeResultEffDose(geant4Result, mcnp6Result, phitsResult)
    writeResultEffDose(geant4Result)
    # writeResultEffDose(mcnp6Result)
    # writeResultEffDose(phitsResult)

    # writeOrganEffectiveDose(geant4Result, mcnp6Result, phitsResult)
    # writeOrganEffectiveDose(geant4Result)
    # writeOrganEffectiveDose(mcnp6Result)
    # writeOrganEffectiveDose(phitsResult)

    # getAnnualDose(geant4Result, mcnp6Result, phitsResult)
