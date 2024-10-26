import pathlib
import shutil
import random

def try_slip_mutation(genome: str) -> str:
    from_idx = random.randint(0, len(genome))
    to_idx = random.randint(0, len(genome) + (from_idx != 0))
    # TODO uncomment one
    #if to_idx < from_idx: # don't allow deletion --- slip insert only
        #return genome
    if to_idx > from_idx: # don't allow insertion --- slip deletion only
        return genome
    return genome[:to_idx] + genome[from_idx:]

def force_slip_mutation(genome: str) -> str:
    result = genome
    if (len(genome) > 100):
        while result == genome or len(result) < 100:
            result = try_slip_mutation(genome)
    return result

lineageFilePath = "detail_MostNumLineageAt50000.dat"

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