import pytest
import os
import pandas as pd
from typing import List
from generateMutationMasks import getMutationMasks

updateToBeAnalyzed = 10000

class Treatment():
    def __init__(self,treatmentPath):
        self.treatmentDir = treatmentPath
        self.runDirectories = []
        self.treatmentName = self.treatmentDir.split('/')[-1]

treatment = Treatment("Baseline-Treatment")
runName = "run_1530"

lineageGenerationIndex = 1

def test_singlePointMutation():
    parent_genome = "wzjagcczyovvccccAvocccpbclqxaycrcccoxurctaiAAlcrcyAccwyccucqcccqcuqylacycctbtcckyccztccctevzvvvvfcaxgab"
    child_genome = "wzjagcczyovvccccAvocccpbclqxaycrcccoxurctaiAAlcrcyAccjyccucqcccqcuqylacycctbtcckyccztccctevzvvvvfcaxgab"

    actualResult = getMutationMasks(
                     parent_genome,
                     child_genome,
                     [53],
                     [],
                     [],
                     treatment,
                     runName,
                     lineageGenerationIndex)
    
    childSourceMap = list(range(len(parent_genome)))
    childSourceMap[53] = -1

    expectedResult = pd.DataFrame(
        {
            "Site": range(len(child_genome)),
            "Treatment": treatment.treatmentName,
            "Run ID": runName,
            "Lineage Generation Index": lineageGenerationIndex,
            "Update Analyzed": updateToBeAnalyzed,
            "CHILD_SOURCE_MAP": childSourceMap,
            "POINT_MUTATION_BOOL_MASK": [i == 53 for i in range(len(child_genome))],
            "SLIP_INSERTION_ORIGIN_BOOL_MASK": [False] * len(parent_genome),
            "SLIP_INSERTION_RESULT_BOOL_MASK": [False] * len(parent_genome),
            "GENOME_CHARACTERS": [child_genome[k] for k in range(len(child_genome))]
        }
    )

    pd.testing.assert_frame_equal(
        actualResult, expectedResult
    )


def test_singleInsertionMutation():
    parent_genome = "wzgcccgab"
    child_genome = "wzgcqcgab"

    actualResult = getMutationMasks(
                     parent_genome,
                     child_genome,
                     [],
                     [],
                     [4],
                     treatment,
                     runName,
                     lineageGenerationIndex)
    

    expectedResult = pd.DataFrame(
        {
            "Site": range(len(child_genome)),
            "Treatment": treatment.treatmentName,
            "Run ID": runName,
            "Lineage Generation Index": lineageGenerationIndex,
            "Update Analyzed": updateToBeAnalyzed,              
            "CHILD_SOURCE_MAP": [0, 1, 2, 3, 4, 5, 6, 7, 8],
            "POINT_MUTATION_BOOL_MASK": [i == 4 for i in range(len(child_genome))],
            "SLIP_INSERTION_ORIGIN_BOOL_MASK": [False] * len(child_genome),
            "SLIP_INSERTION_RESULT_BOOL_MASK": [False] * len(child_genome),                
            "GENOME_CHARACTERS": [child_genome[k] for k in range(len(child_genome))],
        }
    )

    pd.testing.assert_frame_equal(
        actualResult, expectedResult
    )