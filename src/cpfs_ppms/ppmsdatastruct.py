#
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

import numpy as np
from nomad.datamodel.data import (
    ArchiveSection,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
)
from nomad.metainfo import (
    MEnum,
    Quantity,
    Section,
    SubSection,
)


class PPMSData(ArchiveSection):
    """General data section from PPMS"""

    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )

class MPMSData(ArchiveSection):
    '''Maps data section in MPMS'''
    m_def = Section(
        label_quantity='name',
    )
    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )
    map = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')

class MPMSPPMSData(PPMSData):
    '''Data section in MPMS'''
    m_def = Section(
        a_eln=dict(lane_width='600px'),
    )
    time_stamp = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=['*'], description='FILL'
    )
    magnetic_field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    #!TODO: change unit to emu or similar, for now lux=emu
    moment = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    moment_std_err = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    transport_action = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    averaging_time = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    frequency = Quantity(
        type=np.dtype(np.float64), unit='Hz', shape=['*'], description='FILL'
    )
    peak_amplitude = Quantity(
        type=np.dtype(np.float64), unit='mm', shape=['*'], description='FILL'
    )
    center_position = Quantity(
        type=np.dtype(np.float64), unit='mm', shape=['*'], description='FILL'
    )
    lockin_signal_derivative = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    lockin_signal = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    range = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    moment_quad_signal = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    min_temperature = Quantity(
        type=np.dtype(np.float64), unit='K', shape=['*'], description='FILL'
    )
    max_temperature = Quantity(
        type=np.dtype(np.float64), unit='K', shape=['*'], description='FILL'
    )
    min_field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    max_field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    mass = Quantity(
        type=np.dtype(np.float64), unit='g', shape=['*'], description='FILL'
    )
    motor_lag = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    pressure = Quantity(
        type=np.dtype(np.float64), unit='torr', shape=['*'], description='FILL'
    )
    measure_count = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    measurement_number = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    squid_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    motor_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    measure_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    motor_current = Quantity(
        type=np.dtype(np.float64), unit='A', shape=['*'], description='FILL'
    )
    motor_temperature = Quantity(
        type=np.dtype(np.float64), unit='celsius', shape=['*'], description='FILL'
    )
    temperature_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    field_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    chamber_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    chamber_temp = Quantity(
        type=np.dtype(np.float64), unit='K', shape=['*'], description='FILL'
    )
    redirection_state = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    average_temp = Quantity(
        type=np.dtype(np.float64), unit='K', shape=['*'], description='FILL'
    )
    rotation_angle = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    rotator_state = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    dc_moment_fixed_ctr = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    dc_moment_err_fixed_ctr = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    dc_moment_free_ctr = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    dc_moment_err_free_ctr = Quantity(
        type=np.dtype(np.float64), unit='lux', shape=['*'], description='FILL'
    )
    dc_fixed_fit = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    dc_free_fit = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    dc_calculated_center = Quantity(
        type=np.dtype(np.float64), unit='mm', shape=['*'], description='FILL'
    )
    dc_calculated_center_err = Quantity(
        type=np.dtype(np.float64), unit='mm', shape=['*'], description='FILL'
    )
    dc_scan_length = Quantity(
        type=np.dtype(np.float64), unit='mm', shape=['*'], description='FILL'
    )
    dc_scan_time = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    dc_number_of_points = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    dc_squid_drift = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    dc_min_v = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    dc_max_v = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    dc_scans_per_measure = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    sample_position = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    chamber_pressure = Quantity(
        type=np.dtype(np.float64), unit='torr', shape=['*'], description='FILL'
    )
    eto_measurement_mode = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    temperature_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    field_status = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    chamber_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    eto_status_code = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    maps = SubSection(section_def=MPMSData, repeats=True)





class ETOData(ArchiveSection):
    """Data section from Channels in PPMS"""

    m_def = Section(
        label_quantity='name',
    )
    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )
    eto_channel = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')


class ETOChannelData(ArchiveSection):
    """Data section from Channels in PPMS"""

    m_def = Section(
        label_quantity='name',
    )
    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )
    resistance = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    resistance_std_dev = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    phase_angle = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    i_v_current = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    i_v_voltage = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    frequency = Quantity(
        type=np.dtype(np.float64), unit='Hz', shape=['*'], description='FILL'
    )
    averaging_time = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    ac_current = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    dc_current = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    current_ampl = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    in_phase_current = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    quadrature_current = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    gain = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    second_harmonic = Quantity(
        type=np.dtype(np.float64), unit='dB', shape=['*'], description='FILL'
    )
    third_harmonic = Quantity(
        type=np.dtype(np.float64), unit='dB', shape=['*'], description='FILL'
    )


class ETOPPMSData(PPMSData):
    """Data section from PPMS"""

    m_def = Section(
        a_eln=dict(lane_width='600px'),
    )
    time_stamp = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=['*'], description='FILL'
    )
    magnetic_field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    sample_position = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    chamber_pressure = Quantity(
        type=np.dtype(np.float64), unit='torr', shape=['*'], description='FILL'
    )
    eto_measurement_mode = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    temperature_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    field_status = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    chamber_status = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    eto_status_code = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    channels = SubSection(section_def=ETOChannelData, repeats=True)
    eto_channels = SubSection(section_def=ETOData, repeats=True)


class ACTData(ArchiveSection):
    """Data section from Channels in PPMS"""

    m_def = Section(
        label_quantity='name',
    )
    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )
    map = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')


class ACTChannelData(ArchiveSection):
    """Data section from Channels in PPMS"""

    m_def = Section(
        label_quantity='name',
    )
    name = Quantity(
        type=str, description='FILL', a_eln={'component': 'StringEditQuantity'}
    )
    volts = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    v_std_dev = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='ohm/centimeter',
        shape=['*'],
        description='FILL',
    )
    resistivity_std_dev = Quantity(
        type=np.dtype(np.float64),
        unit='ohm/centimeter',
        shape=['*'],
        description='FILL',
    )
    hall = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter**3/coulomb',
        shape=['*'],
        description='FILL',
    )
    hall_std_dev = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter**3/coulomb',
        shape=['*'],
        description='FILL',
    )
    crit_cur = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    crit_cur_std_dev = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    second_harmonic = Quantity(
        type=np.dtype(np.float64), unit='dB', shape=['*'], description='FILL'
    )
    third_harmonic = Quantity(
        type=np.dtype(np.float64), unit='dB', shape=['*'], description='FILL'
    )
    quad_error = Quantity(
        type=np.dtype(np.float64), unit='ohm/cm/rad', shape=['*'], description='FILL'
    )
    drive_signal = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )


class ACTPPMSData(PPMSData):
    """Data section from PPMS"""

    m_def = Section(
        a_eln=dict(lane_width='600px'),
    )
    measurement_type = Quantity(
        type=MEnum(
            'temperature',
            'field',
        ),
        a_eln=ELNAnnotation(
            component='EnumEditQuantity',
        ),
    )
    time_stamp = Quantity(
        type=np.dtype(np.float64), unit='second', shape=['*'], description='FILL'
    )
    status = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin', shape=['*'], description='FILL'
    )
    magnetic_field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    sample_position = Quantity(
        type=np.dtype(np.float64), unit='deg', shape=['*'], description='FILL'
    )
    excitation = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    frequency = Quantity(
        type=np.dtype(np.float64), unit='Hz', shape=['*'], description='FILL'
    )
    act_status = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    act_gain = Quantity(type=np.dtype(np.float64), shape=['*'], description='FILL')
    bridge_1_resistance = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    bridge_2_resistance = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    bridge_3_resistance = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    bridge_4_resistance = Quantity(
        type=np.dtype(np.float64), unit='ohm', shape=['*'], description='FILL'
    )
    bridge_1_excitation = Quantity(
        type=np.dtype(np.float64), unit='microampere', shape=['*'], description='FILL'
    )
    bridge_2_excitation = Quantity(
        type=np.dtype(np.float64), unit='microampere', shape=['*'], description='FILL'
    )
    bridge_3_excitation = Quantity(
        type=np.dtype(np.float64), unit='microampere', shape=['*'], description='FILL'
    )
    bridge_4_excitation = Quantity(
        type=np.dtype(np.float64), unit='microampere', shape=['*'], description='FILL'
    )
    signal_1_vin = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    signal_2_vin = Quantity(
        type=np.dtype(np.float64), unit='V', shape=['*'], description='FILL'
    )
    digital_inputs = Quantity(
        type=np.dtype(np.float64), shape=['*'], description='FILL'
    )
    drive_1_iout = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    drive_2_iout = Quantity(
        type=np.dtype(np.float64), unit='mA', shape=['*'], description='FILL'
    )
    drive_1_ipower = Quantity(
        type=np.dtype(np.float64), unit='watts', shape=['*'], description='FILL'
    )
    drive_2_ipower = Quantity(
        type=np.dtype(np.float64), unit='watts', shape=['*'], description='FILL'
    )
    pressure = Quantity(
        type=np.dtype(np.float64), unit='torr', shape=['*'], description='FILL'
    )
    channels = SubSection(section_def=ACTChannelData, repeats=True)
    maps = SubSection(section_def=ACTData, repeats=True)
