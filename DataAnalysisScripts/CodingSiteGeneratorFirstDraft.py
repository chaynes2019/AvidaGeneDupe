import pandas as pd
import sys
import os

'''
readLineage():
1. Open lineage file
2. Read contents into list of lines
3. Separate into preamble, containing the identifiers for the structure of each line, and the organism 
4. Parse preamble into correctly ordered list of string element identifiers
3. Run parseOrganismInfo() on each line to convert it into a dictionary of organism attributes
4. Return lineage as this list of dictionaries
'''
def readLineage(lineageFile):
    with open(lineageFile, 'r') as f:
        lineageLines = f.readlines()[4:]
        preamble = []
        lineage = []

        for k, line in enumerate(lineageLines):
            if '#' in line:
                preamble.append(line[:-1])
            elif line == '' or line == '\n':
                continue
            else:
                lineage.append(line)

        for k, line in enumerate(preamble):
            itemID = line.split(':')[1]
            preamble[k] = itemID[1:]


        for k, organism in enumerate(lineage):
            if k == 0:
                lineage[k] = parseOrganismInfo(preamble, organism, 1)
            else:
                lineage[k] = parseOrganismInfo(preamble, organism, 0)
    
    return lineage

'''
parseOrganismInfo():
1. Split organism line using spaces
2. Create dictionary by iterating through preamble enumeratively and assigning each preamble element by index
to its corresponding organism string element
3. Return dictionary
'''
def parseOrganismInfo(preamble, organismString, lineageStart):
    preambleCopy = [preamble[k] for k in range(len(preamble))]
    organismStringElements = organismString.split()
    organismAttributes = dict()

    #The first organism in the lineage is lacking parent_muts data; thus, one must readjust
    if lineageStart == 1:
        preambleCopy.pop(2)

    for k, preambleElement in enumerate(preambleCopy):
        organismAttributes[preambleElement] = organismStringElements[k]

    organismAttributes = parseMutationInfo(organismAttributes)

    return organismAttributes

'''
parseMutationInfo():
1. If mutation information is present...
a. Create 3 new lists -- Mutations, Insertions, and Deletions
b. Split mutation information by comma
c. Iterate over each mutation in the split list
d. For each mutation, use the first character to assign the remainder (i.e. the index of the mutation) as an integer to the appropriate list
2. If mutation information is not present...
a. Create 3 new lists -- Mutations, Insertions, and Deletions
b. Leave them empty

The code can directly modify the passed list because Python passes lists by reference
'''
def parseMutationInfo(organismAttributes):
    firstOrganism = False
    
    try:
        mutationInfo = organismAttributes["Mutations from Parent"]
    except(KeyError):
        firstOrganism = True
    
    if not firstOrganism:
        mutations = []
        insertions = []
        deletions = []

        for mutation in mutationInfo.split(','):
            digitCount = [1 for k in range(len(mutation)) if (mutation[k].isdigit())]
            numericalIndexLength = sum(digitCount)

            if mutation[0] == 'M':
                mutations.append(int(mutation[1:(numericalIndexLength + 1)]))
            elif mutation[0] == 'I':
                insertions.append(int(mutation[1:(numericalIndexLength + 1)]))
            elif mutation[0] == 'D':
                deletions.append(int(mutation[1:(numericalIndexLength + 1)]))
            else:
                print("You have something incorrect in the Mutations from Parent object: {}".format(mutation))
        
        organismAttributes["Mutations"] = mutations
        organismAttributes["Insertions"] = insertions
        organismAttributes["Deletions"] = deletions

    else:
        mutations = []
        insertions = []
        deletions = []

        organismAttributes["Mutations"] = mutations
        organismAttributes["Insertions"] = insertions
        organismAttributes["Deletions"] = deletions


'''
findEventPairs():
1. Find Slip event pairs
2. Find Task event pairs
3. Return them as a tuple
'''
def findEventPairs(lineage):
    slipEventPairs = findSlipEventPairs(lineage)
    taskEventPairs = findTaskEventPairs(lineage)

    return (slipEventPairs, taskEventPairs)

'''findSlipEventPairs():
1. Form ordered list of subdictionaries of insertions and deletions pulled from the lineage data
2. For each dictionary (organism), examine the length of the insertions and deletions.
3. If the length of the insertions or deletions > threshold, then add that organism's dictionary from the lineage
to the event pairs list in a list together with the previous organism.
'''


def generateOrganismData():
    return "Beep!"

def writeDataRowToDataFrame():
    return "Keep!"

def generateLineageData():
    stream = os.popen('pwd')
    pwd = stream.read().rstrip()

    runDir = pwd

    updateAtWhichToAnalyze = sys.argv[1]

    lineageDataframe = pd.DataFrame(columns = ["Task Coding Sites",
                                                "Number of Task Coding Sites",
                                                "Number of Unique Coding Sites", 
                                                "Viability Sites", 
                                                "Number of Viability Sites", 
                                                "Genome Length", 
                                                "Fraction Task Coding Sites", 
                                                "Fraction Viability Sites", 
                                                "Ratio of Viability Sites to Coding Sites", 
                                                "Genome",
                                                "Event Type",
                                                "Event Time (Update born of 2nd Organism)"])

    lineageFile = f"detail_MostNumLineageAtUpdate{updateAtWhichToAnalyze}.dat"

    lineage = readLineage(lineageFile)

    (taskEventPairs, duplicationEventPairs) = findEventPairs(lineage)

    for pair in taskEventPairs:
        dataRow = []
        dataRow.append(generateOrganismData(pair[0]))
        dataRow.append(generateOrganismData(pair[1]))
        writeDataRowToDataFrame(dataRow)

    for pair in duplicationEventPairs:
        dataRow = []
        dataRow.append(generateOrganismData(pair[0]))
        dataRow.append(generateOrganismData(pair[1]))
        writeDataRowToDataFrame(dataRow)

    runDirElements = runDir.split('/')
    runName = runDirElements[-1]

    lineageDataframe.to_csv(f"{runName}_LineageDataAtUpdate{updateAtWhichToAnalyze}.csv")

print(readLineage('detail_MostNumLineage.dat'))