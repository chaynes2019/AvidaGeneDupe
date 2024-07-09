import pytest
import os
from enumerateMutants import enumerateMutantInfo

def test_enumeratingAncestor():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert enumerateMutantInfo(testFilePath, "1 -1  0 0 0 0 0 0 0 0 0 1 389 0.249357 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab \n") == ([], [])

def test_allDeletions():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testAllDeletionsFile.dat")
    assert enumerateMutantInfo(testFilePath, "25 10 D25,D40,D60,D61,D62 meep beep sheep \n") == ([], [])

def test_enumeratingBaselineTreatmentPointMutations():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert enumerateMutantInfo(testFilePath, "32862 682 M13cf,M50cb 0 0 0 0 0 0 0 0 0 1 383 0.24282 100 wzcagccccccrcfccccdccccccccccccAccccccccccccccccccbcAcccccccccccccqcccccccccccccccccccccccczvfcaxgab \n") == ([13, 50], [])

def test_enumeratingSlipDuplicateMix():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testSlip-duplicateFile.dat")
    assert enumerateMutantInfo(testFilePath, "25558 598 M79cs,D80c,I91v 0 0 0 0 0 0 0 0 0 1 284 0.341549 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab ") == ([79], [])

def test_enumeratingSlipDuplicateSlipEvent():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testSlip-duplicateFile.dat")
    assert enumerateMutantInfo(testFilePath, "54078 1084 I5c,I6c,I7c,I8c,I9c,I10c,I11c,I12c,I13c,I14c,I15c,I16c,I17c,I18c,I19c,I20c,I21c,I22c,I23c,I24c,I25c,I26c,I27c,I28c,I29c,I30c,I31c,I32c,I33c,I34c,I35o,I36t,I37c,I38c,I39c,I40c,I41c,I42c,I43c,I44c,I45c,I46c,I47c,I48c,I49c,I50c,I51c,I52o,I53c,I54c,I55c,I56c,I57c,I59n,I65y,I79s 1 0 0 0 0 0 0 0 0 1 446 0.343049 156 wzcagccccccccccccccccccccccccccccccotcccccccccccccccoccccccncccccycccccccccccccscccccccccccotcccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n") == ([], [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,59,65,79])