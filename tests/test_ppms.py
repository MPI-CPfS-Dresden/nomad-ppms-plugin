#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import pytest
from nomad.client import normalize_all

# Test ETO functionality

test_files = [
    'tests/data/ETO_Ch1_TMR_Ch2_Hall.dat',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_eto_default(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert (
        parsed_measurement_archive.data.software
        == 'Electrical Transport Option, Release 1.2.0 Build 0'
    )
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 5  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[4].time_stamp) == 3623  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.analyzed_data[0].field) == 1101  # Noqa: PLR2004
    assert parsed_measurement_archive.data.analyzed_data[
        0
    ].carrier_concentration.magnitude == pytest.approx(7.95e27, 1e25)  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 11  # Noqa: PLR2004


# Test ETO LabView functionality

test_files = [
    'tests/data/ETO_LabView_test.xxx',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_eto_labview(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert (
        parsed_measurement_archive.data.software
        == 'Electrical Transport Option - Labview mode'
    )
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 9  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[4].time_stamp) == 37594  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 9  # Noqa: PLR2004


# Test ACT functionality

test_files = [
    'tests/data/ACT_Ch1_Hall_Ch2_TMR.dat',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_act(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert parsed_measurement_archive.data.software == 'ACTRANSPORT,2.0,1.1'
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 8  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[4].time_stamp) == 794  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 8  # Noqa: PLR2004


# Test ACMS funtionality

test_files = [
    'tests/data/ACMS_test.dat',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_acms(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert parsed_measurement_archive.data.software == 'ACMS,1.0,1.1'
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 60  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[0].time_stamp) == 23  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 60  # Noqa: PLR2004


# Test MPMS functionality

test_files = [
    'tests/data/MPMS_test.dat',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_mpms(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert parsed_measurement_archive.data.software == 'MPMS3,1.0,1.1'
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 3  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[0].time_stamp) == 374  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 3  # Noqa: PLR2004


# Test Resistivity functionality

test_files = [
    'tests/data/Resistivity_Ch1_TMR_Ch2_Hall_test.dat',
]
log_levels = ['error', 'critical']


@pytest.mark.parametrize(
    'parsed_measurement_archive, caplog',
    [(file, log_level) for file in test_files for log_level in log_levels],
    indirect=True,
)
def test_normalize_resistivity(parsed_measurement_archive, caplog):
    """
    Tests the normalization of the parsed archive.

    Args:
        parsed_archive (pytest.fixture): Fixture to handle the parsing of archive.
        caplog (pytest.fixture): Fixture to capture errors from the logger.
    """
    normalize_all(parsed_measurement_archive)

    assert parsed_measurement_archive.data.software == 'Resistivity, 2.1, 1.0'
    #  assert len(parsed_measurement_archive.data.steps) == 70 #Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data) == 4  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.data[0].time_stamp) == 199  # Noqa: PLR2004
    assert len(parsed_measurement_archive.data.figures) == 4  # Noqa: PLR2004
