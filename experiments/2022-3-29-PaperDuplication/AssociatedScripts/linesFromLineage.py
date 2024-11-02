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

def do_slip_mutation_deletion(genome: str) -> str:
    num_sites = len(genome)
    if num_sites < 100:
        return genome
    num_available_to_delete = num_sites - 100
    num_to_delete = random.randint(
        int(num_available_to_delete > 0),  # delete at least one, if possible
        num_available_to_delete,
    )
    lo = random.randint(
        0,
        num_sites - num_to_delete,
    )
    hi = lo + num_to_delete
    res = genome[:lo] + genome[hi:]
    assert len(res) <= len(genome)
    assert len(res) >= 100 
    return res

def do_slip_mutation_insertion(genome: str) -> str:
    assert genome
    sites = range(len(genome))
    lo, hi = sorted(random.sample(sites, 2))
    res = genome[:hi] + genome[lo:]
    assert len(res) > len(genome)
    return res


def do_slip_mutation(genome: str) -> str:
    return random.choice(
        [
            do_slip_mutation_deletion,
            do_slip_mutation_insertion,
        ],
    )(genome)

def rewriteLineageFile(runDir):
    lineageFilePath = os.path.join(runDir, f"Timepoint_{desiredUpdateToAnalyze}/data/detail_MostNumLineageAt50000.dat")

    lines = pathlib.Path(lineageFilePath).read_text().splitlines()
    random.seed("".join(lines))

    modifiedLines = []
    for line in lines:
        if "#" not in line and line != "":
            words = line.split(" ")
            genome = words[-2]
            words[-2] = do_slip_mutation(genome)

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