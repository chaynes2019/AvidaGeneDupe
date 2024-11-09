import pathlib
import shutil
import random
import os

desiredUpdateToAnalyze = 50000

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

for subdir in os.listdir(dataDir):
    if subdir not in ["Slip-duplicate"]:
        continue
    treatment = Treatment(os.path.join(dataDir,subdir))
    Treatments.append(treatment)

    for run_dir in os.listdir(treatment.treatmentDir):
        if 'run_' not in run_dir:
            continue
        treatment.runDirectories.append(os.path.join(treatment.treatmentDir, run_dir))

def renameLineageFile(runDir):
    lineageFilePath = os.path.join(runDir, f"Timepoint_{desiredUpdateToAnalyze}/data/detail_MostNumLineageAt50000.dat")
    lineageBackupPath = os.path.join(runDir, f"Timepoint_{desiredUpdateToAnalyze}/data/detail_MostNumLineageAt50000.dat.bak")
    
    lineageFileNewName = os.path.join(runDir, f"Timepoint_{desiredUpdateToAnalyze}/data/detail_MostNumLineageAt50000InsertionBackup.dat")

    os.rename(lineageFilePath, lineageFileNewName)
    os.rename(lineageBackupPath, lineageFilePath)

def renameNewLineageDataFiles(treatmentArray):
    for treatment in treatmentArray:
        treatmentName = treatment.treatmentName
        print(treatmentName)
        
        for runDir in treatment.runDirectories:
            renameLineageFile(runDir)

renameNewLineageDataFiles(Treatments)