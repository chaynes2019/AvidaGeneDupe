import pytest
import os
from enumerateMutants import getID

def test_nullFileIDRetrieval():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testNullFile.dat")
    with pytest.raises(IndexError):
        getID(testFilePath, "")

def test_singleWordIDRetrieval():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testSingleWord.dat")
    assert getID(testFilePath, "apple") == "apple"

def test_baselineTreatmentIDRetrieval():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert getID(testFilePath, "43 104 M11cr 0 0 0 0 0 0 0 0 0 1 388 0.25 100 wzcagccccccrccccccccccccccccccccccccccccccccccccccccAcccccccccccccccccccccccccccccccccccccczvfcaxgab \n") == "43"

def test_SlipDuplicateIDRetrieval():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert getID(testFilePath, "25558 598 M79cs,D80c,I91v 0 0 0 0 0 0 0 0 0 1 284 0.341549 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n") == "25558"
