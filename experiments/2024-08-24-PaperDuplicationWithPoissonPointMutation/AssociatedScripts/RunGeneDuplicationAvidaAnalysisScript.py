from string import Template
import os
import sys

numberOfReplicates = 30
experimentID = sys.argv[1]
updateToAnalyze = sys.argv[2]

validTreatments = ["Slip+", "Slip-", "Slip-_Long"]

treatmentParameters = {"Slip-" : {"Seed Offset" : 1800, "Slip Mutation Probability" : 0.0, "Slip Fill Mode" : 0, "Event File" : "events.cfg"},
                       "Slip+" : {"Seed Offset" : 1830, "Slip Mutation Probability" : 0.05, "Slip Fill Mode" : 0, "Event File" : "events.cfg"},
                       "Slip-_Long" : {"Seed Offset" : 1770, "Slip Mutation Probability" : 0.0, "Slip Fill Mode" : 0, "Event File" : "eventsLongAncestralOrganism.cfg"}}

possibleTreatments = os.listdir(f'/scratch/zamanlh_root/zamanlh0/clhaynes/{experimentID}')

for treatmentInQuestion in possibleTreatments:
    if treatmentInQuestion not in validTreatments:
        continue

    else:
        with open('geneDuplicationDataAnalyzerTemplate.sh', 'r') as templateFile:
            templateString = templateFile.read()
            dataAnalysisScriptTemplate = Template(templateString)

        parameters = treatmentParameters[treatmentInQuestion]

        dataAnalysisScriptString = dataAnalysisScriptTemplate.substitute(treatment=treatmentInQuestion,
                                                                        numReplicates=numberOfReplicates,
                                                                        updateAtWhichToAnalyze=updateToAnalyze,
                                                                        seedOffset=parameters["Seed Offset"],
                                                                        divSlipProb=parameters["Slip Mutation Probability"],
                                                                        slipFillMode=parameters["Slip Fill Mode"],
                                                                        eventFile=parameters["Event File"],
                                                                        experimentalID=experimentID)

        with open(f'geneDuplicationDataAnalyzer_{updateToAnalyze}.sh', 'w') as f:
            print(dataAnalysisScriptString)
            f.write(dataAnalysisScriptString)
        
        os.system(f'bash geneDuplicationDataAnalyzer_{updateToAnalyze}.sh >> spopAnalysisLog.txt')

        os.system(f'rm geneDuplicationDataAnalyzer_{updateToAnalyze}.sh')