import pytest
from testPrep import prepSysPathForTests
prepSysPathForTests()
from CodingSiteGeneratorHPCCCopyFunctions import writeTaskCodingSitesInPandasDataFrame
import CodingSiteGeneratorHPCCCopyFunctions

from unittest.mock import patch, MagicMock
import pandas as pd

CodingSiteGeneratorHPCCCopyFunctions.desiredUpdateToAnalyze = 1000


import os

testDirectory = os.getcwd()


@pytest.fixture
def mock_environment():
    with patch("CodingSiteGeneratorHPCCCopyFunctions.getOrganisms", return_value=["organism1", "organism2"]) as mock_get_organisms, \
         patch("CodingSiteGeneratorHPCCCopyFunctions.getLength", return_value=1000) as mock_get_length, \
         patch("CodingSiteGeneratorHPCCCopyFunctions.getGenome", return_value="ACGTACGTACGTACGT") as mock_get_genome, \
         patch("os.path.join", return_value="dummy_path"):
        yield mock_get_organisms, mock_get_length, mock_get_genome

@pytest.fixture
def mock_treatment():
    treatment = MagicMock()
    treatment.treatmentDataframe = pd.DataFrame(columns=[
        "Run Name", "Task Name", "Timepoint", "Treatment Name", "Coding Sites",
        "Number of Coding Sites", "Number of Unique Coding Sites", "Viability Sites",
        "Number of Viability Sites", "Genome Length", "Fraction of Coding Sites",
        "Fraction of Viability Sites", "Viability to Coding Ratio", "Genome"
    ])
    treatment.treatmentName = "Test Treatment"
    return treatment


def test_write_task_coding_sites_basic(mock_environment, mock_treatment):
    runDir = os.path.join(testDirectory, 'tests/dummy_path')
    taskCodingSites = [[1, 2], [], [3], [], [4, 5], [], [], [6], []]
    viabilitySites = [10, 20, 30]
    numUniqueCodingSites = 5
    writeTaskCodingSitesInPandasDataFrame(mock_treatment, runDir, taskCodingSites, viabilitySites, numUniqueCodingSites)

    runDirElements = runDir.split('/')
    runName = runDirElements[-1]

    assert len(mock_treatment.treatmentDataframe) == 9  # One row per task
    assert mock_treatment.treatmentDataframe.loc["dummy_path,NOT"]["Number of Coding Sites"] == 2
    assert mock_treatment.treatmentDataframe.loc["dummy_path,EQUALS"]["Viability Sites"] == [10, 20, 30]