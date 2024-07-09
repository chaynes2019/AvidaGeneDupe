import pandas as pd
import os
import sys

updateToBeAnalyzed = sys.argv[1]

runDirectories = []
Treatments = []
treatmentParameters = {"Baseline-Treatment":[0.0025, 0.0, 0.0, 0.05, 0.05, 0.0, 0],
"Slip-NOP":[0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 1],
"Slip-duplicate":[0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0],
"Slip-scatter":[0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 5],
"Slip-scramble":[0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 3],
"Slip-random":[0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 2],
"High-Mutation":[0.0025,0.0075,0.0075,0.05,0.05,0.0,0]}
stream = os.popen('pwd')
pwd = stream.read().rstrip()
experimentDir = pwd
dataDir = pwd
experimentName = pwd.split('/')[-1]

class Treatment():
    def __init__(self,treatmentPath):
        self.treatmentDir = treatmentPath
        self.runDirectories = []
        self.treatmentName = self.treatmentDir.split('/')[-1]

dataframe = pd.DataFrame(columns = ["Run ID",
                          "Lineage Generation Index",
                          "Update Analyzed",
                          "Treatment",
                          "Point Mutants",
                          "Slip-Insertion Mutants"])

for subdir in os.listdir(dataDir):
    if subdir not in ['Baseline-Treatment', 'Slip-duplicate']:
        continue
    treatment = Treatment(os.path.join(dataDir,subdir))
    Treatments.append(treatment)

    for run_dir in os.listdir(treatment.treatmentDir):
        if not 'run_' in run_dir:
            continue
        treatment.runDirectories.append(os.path.join(treatment.treatmentDir,run_dir))

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

def enumerateMutantInfo(lineageFile, organism):
    '''Get the ID of the organism to determine
      if mutations should be parsed'''
    orgID = getID(lineageFile, organism)
    
    '''If the organism is the ancestor, it will not have any parent_muts data;
    therefore, define its mutation info as empty and return that value'''
    if(orgID == "1"):
        pointMutants = []
        slipInsertionMutants = []
    
    else:
        pointMutants = []
        slipInsertionMutants = []

        '''Split all the organism line elements into a list.
         Then, split all the individual mutations into a list'''
        mutantInfo = organism.split()[2]
        mutations = mutantInfo.split(',')

        '''Collect all the point mutations and count the number
        of insertion mutations'''
        insertionCounter = 0
        for k in range(len(mutations)):
            mutation = mutations[k]

            digitCount = [1 for j in range(len(mutation)) if (mutation[j].isdigit())]
            numericalIndexLength = sum(digitCount)

            mutationType = mutation[0]
            mutationIndex = mutation[1:(numericalIndexLength + 1)]
            
            try:
                mutationIndex = int(mutationIndex)
            except ValueError:
                raise ValueError(f"The mutation type was {mutationType} and the mutationIndex was {mutationIndex}")

            '''If there is a point mutation, then certainly add
              it to the point mutation collection'''
            if mutationType == 'M':
                pointMutants.append(mutationIndex)

            '''Do not worry about deletion mutations but ignore
            them at present'''
            if mutationType == 'D':
                continue

            '''If there is an insertion mutation, count it.
            The next section of this function distinguishes between
              insertions due to a slip-event and single insertions'''
            if mutationType == 'I':
                insertionCounter += 1
        
        '''If there are more than 10 insertions, it is almost
        certain that a slip-insertion has happened'''
        if insertionCounter > 10:
            for k in range(len(mutations)):
                '''Since the mutations are already sorted by
                index, this should be already sorted in increasing
                order'''
                mutation = mutations[k]

                digitCount = [1 for j in range(len(mutation)) if (mutation[j].isdigit())]
                numericalIndexLength = sum(digitCount)

                mutationType = mutation[0]
                mutationIndex = mutation[1:(numericalIndexLength + 1)]

                try:
                    mutationIndex = int(mutationIndex)
                except ValueError:
                    raise ValueError(f"The mutation type was {mutationType} and the mutationIndex was {mutationIndex}")

                '''Add all insertions to the slip-insertion
                collection -- this will be pruned in the next
                section'''
                if mutationType == 'I':
                    slipInsertionMutants.append(mutationIndex)


    return pointMutants, slipInsertionMutants

def writeMutantsFromLineage(treatment, runName, lineageFile):
    '''1. Read in a lineage data file'''
    organisms = getOrganisms(lineageFile)

    for lineageGenerationIndex in range(len(organisms)):
        '''At each lineage generation index:'''

        '''2. Find point mutants and slip-insertion mutants'''
        pointMutants, slipInsertionMutants = enumerateMutantInfo(lineageFile, organisms[lineageGenerationIndex])

        '''3. Write that in a Pandas dataframe'''
        dataframe.loc[f"{runName}" + f"{lineageGenerationIndex}"] = [runName, lineageGenerationIndex, updateToBeAnalyzed, treatment, pointMutants, slipInsertionMutants]

for treatment in Treatments:
        treatmentName = treatment.treatmentName
        print(treatmentName)
        
        for runDir in treatment.runDirectories:
            lineageFile = os.path.join(runDir, f"Timepoint_{updateToBeAnalyzed}/data/detail_MostNumLineageAt{updateToBeAnalyzed}.dat")
            writeMutantsFromLineage(treatment, runDir, lineageFile)

dataframe.to_csv(f"{experimentDir}/{experimentName}-LineageMutations.csv")