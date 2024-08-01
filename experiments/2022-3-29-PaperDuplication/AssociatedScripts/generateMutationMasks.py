import pandas as pd
from Bio import Align
import itertools as it
from typing import List
import os
import sys

sys.path.append('/home/clhaynes/Documents/AvidaGeneDupe/experiments/2022-3-29-PaperDuplication/AssociatedScripts')

from find_best_alignment import find_best_alignment

def getAlignments():
    pass

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

def getGenome(organism):
    genome = organism.split()[-1]
    return genome

def getSequencePairsFromLineage(lineage_file):
    organisms = getOrganisms(lineage_file)
    
    sequence_pair_collection = []

    for k in range(len(organisms) - 1):
        parent_sequence = getGenome(organisms[k])
        child_sequence = getGenome(organisms[k + 1])

        sequence_pair_collection.append((parent_sequence, child_sequence))
    
    return sequence_pair_collection


def getChildSourceMap(parent_sequence: str,
                      child_sequence: str,
                      alignments: Align.PairwiseAlignments) -> List[int]:
    '''1. Create list of -1 as long as the child_sequence
        This is the starting template'''
    '''2. Collect the alignments tying for the highest score.'''
    '''3. For each alignment in this collection'''
    '''a. Iterate over the subalignments.'''
    '''1. Parse the indices of the child genome involved in the
    subalignment. Parse the indices of the parent genome involved
    in the alignment.'''
    '''2. Using the indices of the child genome subalignment, place
    the indices of the parent genome subalignment into the child
    source map. If the element being replaced is not -1, then raise
    an error.'''

    '''1. Create list of -1 as long as the child_sequence
        This is the starting template'''
    childSourceMap = [-1 for k in range(len(child_sequence))]

    '''2. Collect the alignments tying for the highest score.'''
    #Note: This method presupposes that alignments is already sorted
    #by score in decreasing order
    
    maxScore = 0
    validAlignments = []
    for k in range(len(alignments)):
        '''Define the max score using the first alignment'''
        if k == 0:
            maxScore = alignments[k].score
        
        '''If an alignment is tied with the best alignment,
        then add it for consideration. 
        
        If an alignment is less than the best, then -- 
        because the list of alignments has been sorted
        by score in decreasing order -- end collection'''
        if alignments[k].score == maxScore:
            validAlignments.append(alignments[k])
        elif alignments[k].score < maxScore:
            break
    
        '''3. For each alignment in this collection'''
    for k in range(len(validAlignments)):
        alignmentConsidered = validAlignments[k]

        '''a. Iterate over the subalignments.'''
        #Using len(alignmentConsidered.aligned[0]) because each
        # alignment consists of a 2-tuple; the subtuples of subalignments
        # are what's of interest, and only those subtuples know how long
        #they are
        for idx in range(len(alignmentConsidered.aligned[0])):
            '''Parse the indices of the child genome involved in 
            the subalignment. Parse the indices of the parent genome 
            involved in the alignment.'''
            parentSubAlignment = alignmentConsidered.aligned[0][idx]
            childSubAlignment = alignmentConsidered.aligned[1][idx]

            parent_start = parentSubAlignment[0]
            parent_nonIncEnd = parentSubAlignment[1]

            child_start = childSubAlignment[0]
            child_nonIncEnd = childSubAlignment[1]

            '''Using the indices of the child genome subalignment,
            place the index range of the parent genome subalignment 
            into the child source map. If the element being replaced
            is not -1, then raise an error.'''
            
            for count, j in enumerate(range(child_start, child_nonIncEnd)):
                mappedIndex = range(parent_start, parent_nonIncEnd)[count]

                '''
                if(childSourceMap[j] != -1 and childSourceMap[j] != mappedIndex):
                    raise IndexError(f"getChildSourceMap is trying to overwrite an assigned portion of the map at index {j}: {childSourceMap[j]} to {mappedIndex}")
                else:
                    childSourceMap[j] = mappedIndex
                '''
                
                childSourceMap[j] = mappedIndex

    
    parent_indices = range(len(parent_sequence))

    deletion_mutations = [parent_index for parent_index in parent_indices if parent_index not in childSourceMap]

    slip_insertion_mutations = []
    
    #This complicated series of if-statements is intended to ensure
    #that only bonafide slip-insertions are reported as such, not
    #allowing single insertions of the same instruction type as
    #those they are surrounded by to be reported as slip-insertions
    for k, idx in enumerate(childSourceMap):
        if childSourceMap.count(idx) == 2:
            if idx != -1:
                if k >= 1 and k < len(childSourceMap) - 1:
                    if childSourceMap[k - 1] == idx:
                        childSourceMap[k] = idx
                        childSourceMap[k - 1] = -1
                    elif childSourceMap[k + 1] == idx:
                        childSourceMap[k] = -1
                        childSourceMap[k + 1] = idx
                    else:
                        slip_insertion_mutations.append(k)
                elif k == 0:
                    if childSourceMap[k + 1] == idx:
                        childSourceMap[k] = -1
                        childSourceMap[k + 1] = idx
                    else:
                        slip_insertion_mutations.append(k)
                elif k == len(childSourceMap):
                    if childSourceMap[k - 1] == idx:
                        childSourceMap[k] = idx
                        childSourceMap[k - 1] = -1
                    else:
                        slip_insertion_mutations.append(k)

    insertion_mutations = [k for k in range(len(childSourceMap)) if childSourceMap[k] == -1]

    point_mutations = []
    for child_index, child_value, source_index in zip(it.count(), child_sequence, childSourceMap):
        #Remember: -1 is a valid source_index, but it will direct you
        #to the end of the parent sequence
        if parent_sequence[source_index] != child_value and source_index != -1:
            point_mutations.append(child_index)
            

    return childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations

def getMutationMasks(parent_sequence: str,
                     child_sequence: str,
                     point_mutants: List[int],
                     insertion_mutants: List[int],
                     deletion_mutants: List[int],
                     treatment,
                     runName,
                     lineageGenerationIndex,
                     ) -> pd.DataFrame:
    
    
    childSourceMap = find_best_alignment(parent_sequence, child_sequence)

    pointMutations = []
    for child_index, child_value, source_index in zip(it.count(), child_sequence, childSourceMap):
        #Remember: -1 is a valid source_index, but it will direct you
        #to the end of the parent sequence
        if parent_sequence[source_index] != child_value and source_index != -1:
            pointMutations.append(child_index)
    
    slipInsertionMutations = []
    for k, idx in enumerate(childSourceMap):
        if childSourceMap.count(idx) == 2:
            slipInsertionMutations.append(k)
    
    slipInsertionOriginMutations = []
    slipInsertionResultMutations = []

    if len(slipInsertionMutations) % 2 == 0:
        if len(slipInsertionMutations) > 0:
            #first index of top half = halfwayPoint + 1
            #Therefore, when items are 0-indexed and shifted
            #down by 1, it becomes
            #first index of top half = halfwayPoint
            
            halfwayPoint = len(slipInsertionMutations) // 2
            try:
                slipInsertionOriginMutations = slipInsertionMutations[::halfwayPoint]
                slipInsertionResultMutations = slipInsertionMutations[halfwayPoint::]
            except TypeError as e:
                raise TypeError(f"The slice index was of a wrong type; the index was {halfwayPoint}") from e
    else:
        raise IndexError(f"The method for gathering slip insertion mutations has erred: {slipInsertionMutations}")

    return pd.DataFrame(
        {
            "Site": range(len(child_sequence)),
            "Treatment": treatment.treatmentName,
            "Run ID": runName,
            "Lineage Generation Index": lineageGenerationIndex,
            "Update Analyzed": updateToBeAnalyzed,
            "CHILD_SOURCE_MAP": childSourceMap,
            "POINT_MUTATION_BOOL_MASK": [i in pointMutations for i in range(len(child_sequence))],
            "SLIP_INSERTION_ORIGIN_BOOL_MASK" : [i in slipInsertionOriginMutations for i in range(len(child_sequence))],
            "SLIP_INSERTION_RESULT_BOOL_MASK": [i in slipInsertionResultMutations for i in range(len(child_sequence))],
            "GENOME_CHARACTERS": [child_sequence[k] for k in range(len(child_sequence))]
        }
    )

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

for subdir in os.listdir(dataDir):
    if subdir not in ['Baseline-Treatment', 'Slip-duplicate']:
        continue
    treatment = Treatment(os.path.join(dataDir,subdir))
    Treatments.append(treatment)

    for run_dir in os.listdir(treatment.treatmentDir):
        if not 'run_' in run_dir:
            continue
        treatment.runDirectories.append(os.path.join(treatment.treatmentDir,run_dir))


for treatment in Treatments:
        treatmentName = treatment.treatmentName
        print(treatmentName)
        
        for runDir in treatment.runDirectories:
            lineageFile = os.path.join(runDir, f"Timepoint_{updateToBeAnalyzed}/data/detail_MostNumLineageAt{updateToBeAnalyzed}.dat")

            runDirElements = runDir.split('/')
            runName = runDirElements[-1]

            sequencePairs = getSequencePairsFromLineage(lineageFile)
            print(sequencePairs)

            #Spit out dataframe for the ancestor's origin
            ancestor_sequence = "wzcagcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab"
            mutationMasksDataframe = pd.DataFrame(
                {
                    "Site": range(100),
                    "Treatment": treatment.treatmentName,
                    "Run ID": runName,
                    "Lineage Generation Index": 0,
                    "Update Analyzed": updateToBeAnalyzed,
                    "CHILD_SOURCE_MAP": [-1 for k in range(100)],
                    "POINT_MUTATION_BOOL_MASK": [False] * 100,
                    "SLIP_INSERTION_ORIGIN_BOOL MASK" : [False] * 100,
                    "SLIP_INSERTION_RESULT_BOOL_MASK": [False] * 100,
                    "GENOME_CHARACTERS": [ancestor_sequence[k] for k in range(len(ancestor_sequence))]
                }
            )

            mutationMasksDataframe.to_csv(f"{experimentDir}/{experimentName}-{treatment.treatmentName}-{runName}-Timepoint{updateToBeAnalyzed}-0-LineageMutationMasks.csv")

            for lineageGenerationIdx, sequence_pair in enumerate(sequencePairs):
                lineageGenerationIdx = lineageGenerationIdx + 1

                parentSequence = sequence_pair[0]
                childSequence = sequence_pair[1]

                mutationMasksDataframe = getMutationMasks(parentSequence,
                                                        childSequence,
                                                        [],
                                                        [],
                                                        [],
                                                        treatment,
                                                        runName,
                                                        lineageGenerationIdx
                                                        )

                mutationMasksDataframe.to_csv(f"{experimentDir}/{experimentName}-{treatment.treatmentName}-{runName}-Timepoint{updateToBeAnalyzed}-{lineageGenerationIdx}-LineageMutationMasks.csv")