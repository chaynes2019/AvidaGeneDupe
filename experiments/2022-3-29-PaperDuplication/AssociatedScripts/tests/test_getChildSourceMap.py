import pytest
from Bio import Align
from generateMutationMasks import getChildSourceMap

aligner = Align.PairwiseAligner()
aligner.open_gap_score = -0.25

def test_equalGenomesOfOneLetter():
    parentSequence = "CCCCC"
    childSequence = "CCCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)

    assert childSourceMap == [0, 1, 2, 3, 4]
    assert len(validAlignments) == 1
    assert (parentSubalignment == [0, 5]).all()
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_equalGenomesOfMultipleLetters():
    parentSequence = "ABCDE"
    childSequence = "ABCDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert childSourceMap == [0, 1, 2, 3, 4]
    assert len(validAlignments) == 1
    assert (parentSubalignment == [0, 5]).all()
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])


def test_singlePointMutations():
    parentSequence = "CCCCC"
    childSequence = "CCTCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [0, 5]).all()
    assert childSourceMap == [0, 1, -1, 3, 4]
    assert point_mutations == [2]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_multiplePointMutationsOfSameLetter():
    parentSequence = "CCCCC"
    childSequence = "ACCCA"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [0, 5]).all()
    assert childSourceMap == [-1, 1, 2, 3, -1]
    assert point_mutations == [0, 4]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_manyPointMutations():
    parentSequence = "CCCCC"
    childSequence = "ACACA"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [0, 5]).all()
    assert childSourceMap == [-1, 1, -1, 3, -1]
    assert point_mutations == [0, 2, 4]
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_singleDeletionMutant():
    parentSequence = "CCCCCC"
    childSequence = "CCCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 6
    assert (parentSubalignment == [0, 5]).all()
    assert childSourceMap == [0, 1, 2, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == [5]
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_multipleDeletionsOfSingleLetter():
    parentSequence = "CCCCCC"
    childSequence = "CCCC"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 5
    assert (parentSubalignment == [0, 4]).all()
    assert childSourceMap == [0, 1, 2, 3]
    assert point_mutations == []
    assert deletion_mutations == [4, 5]
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_singleDeletionOfVariedString():
    parentSequence = "ABCDEFG"
    childSequence = "ABCDFG"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [5, 7]).all()
    assert childSourceMap == [0, 1, 2, 3, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == [4]
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_multipleDeletionsOfVariedString():
    parentSequence = "ABCDEFG"
    childSequence = "ACDFG"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [5, 7]).all()
    assert childSourceMap == [0, 2, 3, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == [1, 4]
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_singleInsertionInVariedString():
    parentSequence = "ABCDE"
    childSequence = "ABCFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [3, 5]).all()
    assert childSourceMap == [0, 1, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == [3]
    assert slip_insertion_mutations == set([])

def test_multipleInsertionsInVariedString():
    parentSequence = "ABCDE"
    childSequence = "AXBCFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [3, 5]).all()
    assert childSourceMap == [0, -1, 1, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == [1, 4]
    assert slip_insertion_mutations == set([])

def test_mixedInsertionAndDeletionInVariedString():
    parentSequence = "ABCDE"
    childSequence = "ACFDE"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [3, 5]).all()
    assert childSourceMap == [0, 2, -1, 3, 4]
    assert point_mutations == []
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == set([])

def test_mixedInsertionDeletionPointMutationsInVariedString():
    parentSequence = "ABCDEF"
    childSequence = "ACXDET"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [3, 6]).all()
    assert childSourceMap == [0, 2, -1, 3, 4, -1]
    assert point_mutations == [5]
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == set([])

def test_pointMutationMatchesEnd():
    parentSequence = "ABCDEF"
    childSequence = "ACFDEF"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [3, 6]).all()
    assert childSourceMap == [0, 2, -1, 3, 4, 5]
    assert point_mutations == []
    assert deletion_mutations == [1]
    assert insertion_mutations == [2]
    assert slip_insertion_mutations == set([])

def test_realBaseline_TreatmentExample():
    parentSequence = "wzcakcuppstrrtkccisdxwycywmcbccpqgcccctpcibecqeccccscqcycuyoccpccecnuccccecsycccqmyAdcluccczvvfvvxgab"
    childSequence = "wzcakcuppstrrtkccisdxwycywmcbccpqgcccctpcibecqeccccscqcycuyoccpccecnuccccecsycccqmyAdcluccczvvvvxgab"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 1
    assert (parentSubalignment == [95, 101]).all()
    assert childSourceMap == [k for k in range(94)] + [k + 95 for k in range(6)]
    assert point_mutations == []
    assert deletion_mutations == [94]
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([])

def test_insertionNextToOneOfSameInstruction():
    parentSequence = "AFD"
    childSequence = "AFFD"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 2
    assert (parentSubalignment == [2, 3]).all()
    assert childSourceMap == [0, 1, 1, 2]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([1])

def test_slipDuplication():
    parentSequence = "ABDEFWX"
    childSequence = "ABDEFDEFWX"
    alignmentList = aligner.align(parentSequence, childSequence)

    childSourceMap, validAlignments, parentSubalignment, point_mutations, deletion_mutations, insertion_mutations, slip_insertion_mutations = getChildSourceMap(parentSequence, 
                                                                                                                                      childSequence, 
                                                                                                                                      alignmentList)
    
    assert len(validAlignments) == 4
    assert (parentSubalignment == [5, 7]).all()
    assert childSourceMap == [0, 1, 2, 3, 4, 2, 3, 4, 5, 6]
    assert point_mutations == []
    assert deletion_mutations == []
    assert insertion_mutations == []
    assert slip_insertion_mutations == set([2, 3, 4])