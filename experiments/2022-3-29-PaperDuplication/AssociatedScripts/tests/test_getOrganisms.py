import pytest
import os
from enumerateMutants import getOrganisms

def test_getSingleWord():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testSingleWord.dat")
    assert getOrganisms(testFilePath) == ["apple"]

def test_getNullFile():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testNullFile.dat")
    with pytest.raises(ValueError):
        getOrganisms(testFilePath)

def test_normalLineageFile():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testNormalLineage.dat")
    assert getOrganisms(testFilePath) == ["apple"]
