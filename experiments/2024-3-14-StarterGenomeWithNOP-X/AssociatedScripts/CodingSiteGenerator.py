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
            lineage[k] = parseOrganismInfo(preamble, organism)
    
    return lineage

'''
parseOrganismInfo():
1. Split organism line using spaces
2. Create dictionary by iterating through preamble enumeratively and assigning each preamble element by index
to its corresponding organism string element
3. Return dictionary
'''
def parseOrganismInfo(preamble, organismString):
    organismStringElements = organismString.split()
    organismAttributes = dict()

    for k, preambleElement in enumerate(preamble):
        organismAttributes[preambleElement] = organismStringElements[k]

    return organismAttributes

def findEventPairs():
    return "Sneap!"

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

print(readLineage('detail-MostNumLineageAt10000.dat'))