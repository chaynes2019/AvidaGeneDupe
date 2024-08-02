import pytest
from Bio import Align
from generateMutationMasks import getChildSourceMap

aligner = Align.PairwiseAligner()
aligner.open_gap_score = -0.25

def test_equalGenomesOfOneLetter():
    parentSequence = "CCCCC"
    childSequence = "CCCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)

    assert childSourceMap == [0, 1, 2, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_equalGenomesOfMultipleLetters():
    parentSequence = "ABCDE"
    childSequence = "ABCDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == []
    assert childSourceMap == [0, 1, 2, 3, 4]


def test_singlePointMutations():
    parentSequence = "CCCCC"
    childSequence = "CCTCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, -1, 3, 4]
    assert point_mutations == [2]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_multiplePointMutationsOfSameLetter():
    parentSequence = "CCCCC"
    childSequence = "ACCCA"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [-1, 1, 2, 3, -1]
    assert point_mutations == [0, 4]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_manyPointMutations():
    parentSequence = "CCCCC"
    childSequence = "ACACA"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [-1, 1, -1, 3, -1]
    assert point_mutations == [0, 2, 4]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_singleDeletionMutant():
    parentSequence = "CCCCCC"
    childSequence = "CCCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == [5]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_multipleDeletionsOfSingleLetter():
    parentSequence = "CCCCCC"
    childSequence = "CCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3]
    assert point_mutations == []
    assert deletion_mutations == [4, 5]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_singleDeletionOfVariedString():
    parentSequence = "ABCDEFG"
    childSequence = "ABCDFG"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    assert childSourceMap == [0, 1, 2, 3, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == [4]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_multipleDeletionsOfVariedString():
    parentSequence = "ABCDEFG"
    childSequence = "ACDFG"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 2, 3, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == [1, 4]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_singleInsertionInVariedString():
    parentSequence = "ABCDE"
    childSequence = "ABCFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    
    assert childSourceMap == [0, 1, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == [3]
    assert slip_insertion_mutations == []

def test_multipleInsertionsInVariedString():
    parentSequence = "ABCDE"
    childSequence = "AXBCFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, -1, 1, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == [1, 4]
    assert slip_insertion_mutations == []

def test_mixedInsertionAndDeletionInVariedString():
    parentSequence = "ABCDE"
    childSequence = "ACFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    
    assert childSourceMap == [0, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == []

def test_mixedInsertionDeletionPointMutationsInVariedString():
    parentSequence = "ABCDEF"
    childSequence = "ACXDET"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 2, -1, 3, 4, -1]
    assert point_mutations == [5]
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == []

def test_pointMutationMatchesEnd():
    parentSequence = "ABCDEF"
    childSequence = "ACFDEF"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 2, -1, 3, 4, 5]
    assert point_mutations == []
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == []

def test_realBaseline_TreatmentExample():
    parentSequence = "wzcakcuppstrrtkccisdxwycywmcbccpqgcccctpcibecqeccccscqcycuyoccpccecnuccccecsycccqmyAdcluccczvvfvvxgab"
    childSequence = "wzcakcuppstrrtkccisdxwycywmcbccpqgcccctpcibecqeccccscqcycuyoccpccecnuccccecsycccqmyAdcluccczvvvvxgab"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [k for k in range(94)] + [k + 95 for k in range(6)]
    assert point_mutations == []
    assert deletion_mutations == [94]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_insertionNextToOneOfSameInstruction():
    parentSequence = "AFD"
    childSequence = "AFFD"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, -1, 1, 2]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == [1]
    assert slip_insertion_mutations == []

def test_slipDuplication():
    parentSequence = "ABDEFWX"
    childSequence = "ABDEFDEFWX"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 4, 2, 3, 4, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == [2, 3, 4, 5, 6, 7]

def test_slipDeletion():
    parentSequence = "ABDEFYUPWX"
    childSequence = "ABDEWX"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 8, 9]
    assert point_mutations == []
    assert deletion_mutations == [4, 5, 6, 7]
    assert insertion_mutations == []
    assert slip_insertion_mutations == []

def test_aaaaaStrangePhenomena():
    parentSequence = "AAAAWXYZ"
    childSequence = "AAAAAAAAWXYZ"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == [0, 1, 2, 3, 4, 5, 6, 7]

def test_criscrossingRun_1609_OneSite_ManyChildren():
    parentSequence = "wzcagcccccccccccccccccccccccccccccocccccccccccccclcccccccccccccccccccccccccccdcccccccccccczvvfcaxgab"
    childSequence = "wzcagcccccccccccccccccccccccccccccocccccccccccccclcccccccccccccccccccccccccccdcccscccccccczvvfcaxgab"

    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == [0, 1, 2, 3, 4, 5, 6, 7]

def test_spacedDuplication():
    parentSequence = "wzcagcclacccneccccccccncqccccccccutcccvccccccccccccefcccccccacccycccccccccncccccccccccqccccccmcucccccccccccccnccceccccccccacccycccccccacccycccccccccncccccccccccqccccccmcuccccccccccccccccceccccccccacccycccccccccccccccruccqccccccccycccccccccncccccccccccqccccccmcucccccccuccccccccceccccccccacccycccccccccccccccruccqccccccczvfcaxgab"
    childSequence =  "wzcagcclacccneccccccccncqccccccccutcccvcccccccccckcefcccccccacccycccccccccncccccccccccqccccccmcucccccccccccccnccceccccccccacccycccccccacccycccccccccncccccccccccqccccccmcuccccccccccccccccceccccccccacccyccccccccccccccruccqccccccccycccccccccncccccccccccqccccycccccccccccccccruccqccccccccycccccccccncccccccccccqccccccmcucccccccuccccccccceccccccccacccycccccccccccccccruccqccccccczvfcaxgab"

    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    expectedChildSourceMap = [k for k in range(215)] + [-1 for k in range(37)] + [197 + k for k in range(131)]
    expectedChildSourceMap[49] = -1


    assert childSourceMap == expectedChildSourceMap
    assert point_mutations == [49]
    assert deletion_mutations == []
    assert insertion_mutations == [215 + k for k in range(37)]
    assert slip_insertion_mutations == [197 + k for k in range(18)] + [252 + k for k in range(18)]