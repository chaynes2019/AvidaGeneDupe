import pandas as pd
import sys
import os

def readLineage(lineageFile):
    with open(lineageFile, 'r') as f:
        lineageLines = f.readlines()

        for k, organism in enumerate(lineageLines):
            lineageLines[k] = parseOrganismInfo(organism)
    
    lineage = lineageLines

    return lineage

def parseOrganismInfo(organismString):
    return "Meep!"

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

    lineageDataframe.to_csv(f"{treatmentName}_{runName}_LineageDataAtUpdate{updateAtWhichToAnalyze}.csv")

generateLineageData()