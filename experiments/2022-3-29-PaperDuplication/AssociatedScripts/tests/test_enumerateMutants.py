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
    assert enumerateMutantInfo(testFilePath, "25 D25,D40,D60,D61,D62 meep beep sheep \n") == ([], [])

def test_enumeratingBaselineTreatmentMix():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert enumerateMutantInfo(testFilePath, )