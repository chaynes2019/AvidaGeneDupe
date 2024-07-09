import pandas as pd
import sys

updateToBeAnalyzed = sys.argv[1]

dataframe = pd.DataFrame(["Run ID",
                          "Lineage Generation Index",
                          "Update Analyzed",
                          "Treatment"
                          "Point Mutants",
                          "Slip-Insertion Mutants"])

def getOrganisms(lineageFile):
    with open(lineageFile,'r') as datFile:
        lines = datFile.readlines()
        initialOrgPos = 0
        for k,line in enumerate(lines):
            if (line[0] != '') & (line[0] != '#') & (line[0] != '\n'):
                initialOrgPos = k
                break
            else:
                continue

        organisms = []
        for i in range(initialOrgPos,len(lines)):
            if(lines[i] != ''):
                organisms.append(lines[i])
            else:
                continue

        if organisms == []:
            raise ValueError(f"The lineage file has no organisms; check file provided: {lineageFile}")
        else:
            return organisms

def getID(lineageFile, organism):
    '''1. Split organism line by spaces'''
    organismElements = organism.split()

    '''2. The ID is the first element in the organism line'''
    try:
        id = organismElements[0]
    except(IndexError) as e:
        raise IndexError(f"There are no characters in the organism line in {lineageFile}") from e

    return id

def enumerateMutants(lineageFile, organism):
    pass

def writeMutantsFromLineage(treatment, runName, lineageFile):
    '''1. Read in a lineage data file'''
    organisms = getOrganisms(lineageFile)

    for lineageGenerationIndex in range(len(organisms)):
        '''At each lineage generation index:'''

        '''2. Find point mutants and slip-insertion mutants'''
        pointMutants, slipInsertionMutants = enumerateMutants(organisms[lineageGenerationIndex])

        '''3. Write that in a Pandas dataframe'''
        dataframe[lineageGenerationIndex] = [runName, lineageGenerationIndex, updateToBeAnalyzed, treatment, pointMutants, slipInsertionMutants]
