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

def test_baselineTreatmentLineageFile():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testBaseline-TreatmentFile.dat")
    assert getOrganisms(testFilePath) == ["1 -1  0 0 0 0 0 0 0 0 0 1 389 0.249357 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "18 78 M52cA 0 0 0 0 0 0 0 0 0 1 389 0.249357 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccAcccccccccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "43 104 M11cr 0 0 0 0 0 0 0 0 0 1 388 0.25 100 wzcagccccccrccccccccccccccccccccccccccccccccccccccccAcccccccccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "783 160 M18cd 0 0 0 0 0 0 0 0 0 1 387 0.250646 100 wzcagccccccrccccccdcccccccccccccccccccccccccccccccccAcccccccccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "6091 253 M66cq 0 0 0 0 0 0 0 0 0 1 387 0.250646 100 wzcagccccccrccccccdcccccccccccccccccccccccccccccccccAcccccccccccccqcccccccccccccccccccccccczvfcaxgab \n",
                                          "28461 614 M31cA 0 0 0 0 0 0 0 0 0 1 387 0.250646 100 wzcagccccccrccccccdccccccccccccAccccccccccccccccccccAcccccccccccccqcccccccccccccccccccccccczvfcaxgab \n",
                                          "32862 682 M13cf,M50cb 0 0 0 0 0 0 0 0 0 1 383 0.24282 100 wzcagccccccrcfccccdccccccccccccAccccccccccccccccccbcAcccccccccccccqcccccccccccccccccccccccczvfcaxgab \n",
                                          "47825 924 M30cu 0 0 0 0 0 0 0 0 0 1 383 0.24282 100 wzcagccccccrcfccccdcccccccccccuAccccccccccccccccccbcAcccccccccccccqcccccccccccccccccccccccczvfcaxgab \n"]

def test_slipDuplicateLineageFile():
    wd = os.getcwd()
    testFilePath = os.path.join(wd, "tests/TestFiles/testSlip-duplicateFile.dat")
    assert getOrganisms(testFilePath) == ["1 -1  0 0 0 0 0 0 0 0 0 1 389 0.249357 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "4911 237 M59cn 0 0 0 0 0 0 0 0 0 1 388 0.25 100 wzcagccccccccccccccccccccccccccccccccccccccccccccccccccccccnccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "9148 315 M52co 0 0 0 0 0 0 0 0 0 1 387 0.250646 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccoccccccnccccccccccccccccccccccccccccccczvfcaxgab \n",
                                          "9969 329 M65cy 0 0 0 0 0 0 0 0 0 1 386 0.251295 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccoccccccncccccyccccccccccccccccccccccccczvfcaxgab \n",
                                          "25558 598 M79cs,D80c,I91v 0 0 0 0 0 0 0 0 0 1 284 0.341549 100 wzcagcccccccccccccccccccccccccccccccccccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n",
                                          "33066 728 M36ct 0 0 0 0 0 0 0 0 0 1 283 0.342756 100 wzcagccccccccccccccccccccccccccccccctcccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n",
                                          "38532 821 M35cl 0 0 0 0 0 0 0 0 0 1 283 0.342756 100 wzcagccccccccccccccccccccccccccccccltcccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n",
                                          "47441 973 M35lo 0 0 0 0 0 0 0 0 0 1 283 0.342756 100 wzcagccccccccccccccccccccccccccccccotcccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n",
                                          "54078 1084 I5c,I6c,I7c,I8c,I9c,I10c,I11c,I12c,I13c,I14c,I15c,I16c,I17c,I18c,I19c,I20c,I21c,I22c,I23c,I24c,I25c,I26c,I27c,I28c,I29c,I30c,I31c,I32c,I33c,I34c,I35o,I36t,I37c,I38c,I39c,I40c,I41c,I42c,I43c,I44c,I45c,I46c,I47c,I48c,I49c,I50c,I51c,I52o,I53c,I54c,I55c,I56c,I57c,I59n,I65y,I79s 1 0 0 0 0 0 0 0 0 1 446 0.343049 156 wzcagccccccccccccccccccccccccccccccotcccccccccccccccoccccccncccccycccccccccccccscccccccccccotcccccccccccccccoccccccncccccycccccccccccccscccccccccczvvfcaxgab \n"]
