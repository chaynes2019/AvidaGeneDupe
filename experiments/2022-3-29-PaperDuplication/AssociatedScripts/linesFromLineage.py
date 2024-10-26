import pathlib
import shutil
import random
import os
import numpy as np
import pandas as pd
import sys
import uuid

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

def try_slip_mutation(genome: str) -> str:
    from_idx = random.randint(0, len(genome))
    to_idx = random.randint(0, len(genome) + (from_idx != 0))
    # TODO uncomment one
    if to_idx < from_idx: # don't allow deletion --- slip insert only
        return genome
    #if to_idx > from_idx: # don't allow insertion --- slip deletion only
        #return genome
    return genome[:to_idx] + genome[from_idx:]

def force_slip_mutation(genome: str) -> str:
    result = genome
    if (len(genome) > 100):
        while result == genome or len(result) < 100:
            result = try_slip_mutation(genome)
    return result

def rewriteLineageFile(runDir):
    lineageFilePath = os.path.join(runDir, f"Timepoint_{desiredUpdateToAnalyze}/data/detail_MostNumLineageAt50000.dat")

    lines = pathlib.Path(lineageFilePath).read_text().splitlines()
    random.seed("".join(lines))

    modifiedLines = []
    for line in lines:
        if "#" not in line and line != "":
            words = line.split(" ")
            genome = words[-2]
            words[-2] = force_slip_mutation(genome)

            line = " ".join(words)

        
        modifiedLines.append(line)

    shutil.copy(lineageFilePath, lineageFilePath + ".bak")

    with open(lineageFilePath, 'w') as f:
        for line in modifiedLines:
            f.write(line + "\n")

def writeNewLineageDataFiles(treatmentArray):
    for treatment in treatmentArray:
        treatmentName = treatment.treatmentName
        print(treatmentName)
        
        for runDir in treatment.runDirectories:
            rewriteLineageFile(runDir)

writeNewLineageDataFiles(Treatments)