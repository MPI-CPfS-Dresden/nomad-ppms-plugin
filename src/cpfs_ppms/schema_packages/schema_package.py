# ruff: noqa: E501

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

import re
from io import StringIO
from typing import (
    TYPE_CHECKING,
)

import pandas as pd
from nomad.datamodel.data import (
    EntryData,
)
from nomad.datamodel.metainfo.plot import PlotSection
from nomad.search import search
from nomad_ppms_plugin.cpfsppmsdatastruct import (
    CPFSPPMSACMSMeasurement,
    CPFSPPMSACTMeasurement,
    CPFSPPMSETOMeasurement,
    CPFSPPMSMPMSMeasurement,
    CPFSSample,
)
from nomad_ppms_plugin.ppmsfunctions import (
    get_fileopentime,
    get_ppms_steps_from_data,
    split_ppms_data_acms,
    split_ppms_data_act,
    split_ppms_data_eto,
    split_ppms_data_mpms,
)
from structlog.stdlib import (
    BoundLogger,
)

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.config import config
from nomad.metainfo import SchemaPackage

configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.schema_packages:schema_entry_point_eto_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.schema_packages:schema_entry_point_eto_labview'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.schema_packages:schema_entry_point_act_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.schema_packages:schema_entry_point_mpms_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.schema_packages:schema_entry_point_acms_default'
)

m_package_ppms_eto_default = SchemaPackage()


class CPFSPPMSETOMeasurementDefault(CPFSPPMSETOMeasurement, PlotSection, EntryData):
    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            # For automatic step discovery, some parameters are needed:
            if not self.temperature_tolerance:
                self.temperature_tolerance = 0.2
            if not self.field_tolerance:
                self.field_tolerance = 5.0

            with archive.m_context.raw_file(self.data_file, 'r') as file:
                data = file.read()

            header_match = re.search(r'\[Header\](.*?)\[Data\]', data, re.DOTALL)
            header_section = header_match.group(1).strip()
            header_lines = header_section.split('\n')

            while self.samples:
                self.m_remove_sub_section(CPFSPPMSETOMeasurementDefault.samples, 0)

            for i in ['1', '2']:
                sample_headers = [
                    line
                    for line in header_lines
                    if line.startswith('INFO') and 'SAMPLE' + i + '_' in line
                ]
                if sample_headers:
                    sample = CPFSSample()
                    for line in sample_headers:
                        parts = re.split(r',\s*', line)
                        key = parts[-1].lower().replace('sample' + i + '_', '')
                        if key == 'material':
                            for line2 in parts[1:-1]:
                                if line2.startswith('l='):
                                    setattr(
                                        sample,
                                        'length',
                                        float(line2.strip('l=').strip('mm')) / 1000.0,
                                    )
                                if line2.startswith('w='):
                                    setattr(
                                        sample,
                                        'width',
                                        float(line2.strip('w=').strip('mm')) / 1000.0,
                                    )
                                if line2.startswith('t='):
                                    setattr(
                                        sample,
                                        'depth',
                                        float(line2.strip('t=').strip('mm')) / 1000.0,
                                    )
                        if key == 'comment':
                            setattr(sample, key, ', '.join(parts[1:-1]))
                            # ids="_".join(parts[1].split("_")[1:3])
                            ids = parts[1].split('_')[1]
                            logger.info(ids)
                            search_result = search(
                                owner='user',
                                query={
                                    'results.eln.sections:any': ['CPFSCrystal'],
                                    'results.eln.names:any': [ids + r'*'],
                                },
                                user_id=archive.metadata.main_author.user_id,
                            )
                            if len(search_result.data) > 0:
                                sample.sample_id = f'../uploads/{search_result.data[0]["upload_id"]}/archive/{search_result.data[0]["entry_id"]}#data'
                                sample.name = search_result.data[0][
                                    'search_quantities'
                                ][0]['str_value']
                            else:
                                logger.warning(
                                    "The sample given in the header could not be found and couldn't be referenced."
                                )
                        elif hasattr(sample, key):
                            setattr(sample, key, ', '.join(parts[1:-1]))
                    if not sample.length:
                        logger.info(
                            'Length for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'length', 1.0)
                    if not sample.width:
                        logger.info(
                            'Width for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'width', 1.0)
                    if not sample.depth:
                        logger.info(
                            'Depth for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'depth', 1.0)
                    self.m_add_sub_section(
                        CPFSPPMSETOMeasurementDefault.samples, sample
                    )

            for line in header_lines:
                if line.startswith('FILEOPENTIME'):
                    if hasattr(self, 'datetime'):
                        iso_date = get_fileopentime(line)
                        if iso_date == 'Not found.':
                            logger.error(
                                'FILEOPENTIME not understood. Check the format.'
                            )
                        else:
                            setattr(self, 'datetime', iso_date)
                if line.startswith('BYAPP'):
                    if hasattr(self, 'software'):
                        setattr(self, 'software', line.replace('BYAPP,', '').strip())
                if line.startswith('TEMPERATURETOLERANCE'):
                    if hasattr(self, 'temperature_tolerance'):
                        setattr(
                            self,
                            'temperature_tolerance',
                            float(line.replace('TEMPERATURETOLERANCE,', '').strip()),
                        )
                if line.startswith('FIELDTOLERANCE'):
                    if hasattr(self, 'field_tolerance'):
                        setattr(
                            self,
                            'field_tolerance',
                            float(line.replace('FIELDTOLERANCE,', '').strip()),
                        )

            data_section = header_match.string[header_match.end() :]
            data_section = data_section.replace(',Field', ',Magnetic Field')
            data_buffer = StringIO(data_section)
            data_df = pd.read_csv(
                data_buffer,
                header=0,
                skipinitialspace=True,
                sep=',',
                engine='python',
            )

            all_steps, runs_list = get_ppms_steps_from_data(
                data_df, self.temperature_tolerance, self.field_tolerance
            )

            if not self.sequence_file:
                self.steps = all_steps

            logger.info('Parsing ETO measurement.')
            logger.info(runs_list)
            self.data = split_ppms_data_eto(data_df, runs_list)

        super().normalize(archive, logger)


m_package_ppms_eto_default.__init_metainfo__()


m_package_ppms_eto_labview = SchemaPackage()


class CPFSPPMSETOMeasurementLabview(CPFSPPMSETOMeasurement, PlotSection, EntryData):
    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        self.software = 'Electrical Transport Option - Labview mode'

        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            # For automatic step discovery, some parameters are needed:
            if not self.temperature_tolerance:
                self.temperature_tolerance = 0.05
            if not self.field_tolerance:
                self.field_tolerance = 5.0

            with archive.m_context.raw_file(self.data_file, 'r') as file:
                data_section = file.read()

            data_section = re.sub('\t', ', ', data_section)
            data_section = data_section.replace('time_passed', 'time stamp')
            data_section = data_section.replace(
                'mercury_IPS:_Field', 'Magnetic Field (Oe)'
            )
            data_section = data_section.replace(
                'SR830_#2:_X-value', 'Resistance Ch1 (Ohms)'
            )
            data_section = data_section.replace(
                'SR830_#3:_X-value', 'Resistance Ch2 (Ohms)'
            )
            data_section = data_section.replace(
                'mercury_temp_ctrl:_VTI:_Temperature', 'Temperature (K)'
            )

            data_buffer = StringIO(data_section)
            data_df = pd.read_csv(
                data_buffer,
                header=0,
                skipinitialspace=True,
                sep=',',
                comment=';',
                engine='python',
            )

            data_df['Magnetic Field (Oe)'] = data_df['Magnetic Field (Oe)'] * 10000

            all_steps, runs_list = get_ppms_steps_from_data(
                data_df, self.temperature_tolerance, self.field_tolerance
            )

            if not self.sequence_file:
                self.steps = all_steps

            self.data = split_ppms_data_eto(data_df, runs_list)

        super().normalize(archive, logger)


m_package_ppms_eto_labview.__init_metainfo__()

m_package_ppms_act_default = SchemaPackage()


class CPFSPPMSACTMeasurementDefault(CPFSPPMSACTMeasurement, PlotSection, EntryData):
    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            # For automatic step discovery, some parameters are needed:
            if not self.temperature_tolerance:
                self.temperature_tolerance = 0.2
            if not self.field_tolerance:
                self.field_tolerance = 5.0

            with archive.m_context.raw_file(self.data_file, 'r') as file:
                data = file.read()

            header_match = re.search(r'\[Header\](.*?)\[Data\]', data, re.DOTALL)
            header_section = header_match.group(1).strip()
            header_lines = header_section.split('\n')

            while self.samples:
                self.m_remove_sub_section(CPFSPPMSACTMeasurementDefault.samples, 0)

            for i in ['1', '2']:
                sample_headers = [
                    line
                    for line in header_lines
                    if line.startswith('INFO') and 'SAMPLE' + i + '_' in line
                ]
                if sample_headers:
                    sample = CPFSSample()
                    for line in sample_headers:
                        parts = re.split(r',\s*', line)
                        key = parts[-1].lower().replace('sample' + i + '_', '')
                        if key == 'material':
                            for line2 in parts[1:-1]:
                                if line2.startswith('l='):
                                    setattr(
                                        sample,
                                        'length',
                                        float(line2.strip('l=').strip('mm')) / 1000.0,
                                    )
                                if line2.startswith('w='):
                                    setattr(
                                        sample,
                                        'width',
                                        float(line2.strip('w=').strip('mm')) / 1000.0,
                                    )
                                if line2.startswith('t='):
                                    setattr(
                                        sample,
                                        'depth',
                                        float(line2.strip('t=').strip('mm')) / 1000.0,
                                    )
                        if key == 'comment':
                            setattr(sample, key, ', '.join(parts[1:-1]))
                            # ids="_".join(parts[1].split("_")[1:3])
                            ids = parts[1].split('_')[1]
                            logger.info(ids)
                            search_result = search(
                                owner='user',
                                query={
                                    'results.eln.sections:any': ['CPFSCrystal'],
                                    'results.eln.names:any': [ids + r'*'],
                                },
                                user_id=archive.metadata.main_author.user_id,
                            )
                            if len(search_result.data) > 0:
                                sample.sample_id = f'../uploads/{search_result.data[0]["upload_id"]}/archive/{search_result.data[0]["entry_id"]}#data'
                                sample.name = search_result.data[0][
                                    'search_quantities'
                                ][0]['str_value']
                            else:
                                logger.warning(
                                    "The sample given in the header could not be found and couldn't be referenced."
                                )
                        elif hasattr(sample, key):
                            setattr(sample, key, ', '.join(parts[1:-1]))
                    if not sample.length:
                        logger.info(
                            'Length for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'length', 1.0)
                    if not sample.width:
                        logger.info(
                            'Width for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'width', 1.0)
                    if not sample.depth:
                        logger.info(
                            'Depth for sample '
                            + str(i)
                            + ' not found, set to 1. Value of resistivity will be wrong.'
                        )
                        setattr(sample, 'depth', 1.0)
                    self.m_add_sub_section(
                        CPFSPPMSACTMeasurementDefault.samples, sample
                    )

            for line in header_lines:
                if line.startswith('FILEOPENTIME'):
                    if hasattr(self, 'datetime'):
                        iso_date = get_fileopentime(line)
                        if iso_date == 'Not found.':
                            logger.error(
                                'FILEOPENTIME not understood. Check the format.'
                            )
                        else:
                            setattr(self, 'datetime', iso_date)
                if line.startswith('BYAPP'):
                    if hasattr(self, 'software'):
                        setattr(self, 'software', line.replace('BYAPP,', '').strip())
                if line.startswith('TEMPERATURETOLERANCE'):
                    if hasattr(self, 'temperature_tolerance'):
                        setattr(
                            self,
                            'temperature_tolerance',
                            float(line.replace('TEMPERATURETOLERANCE,', '').strip()),
                        )
                if line.startswith('FIELDTOLERANCE'):
                    if hasattr(self, 'field_tolerance'):
                        setattr(
                            self,
                            'field_tolerance',
                            float(line.replace('FIELDTOLERANCE,', '').strip()),
                        )

            data_section = header_match.string[header_match.end() :]
            data_section = data_section.replace(',Field', ',Magnetic Field')
            data_buffer = StringIO(data_section)
            data_df = pd.read_csv(
                data_buffer,
                header=0,
                skipinitialspace=True,
                sep=',',
                engine='python',
            )

            all_steps, runs_list = get_ppms_steps_from_data(
                data_df, self.temperature_tolerance, self.field_tolerance
            )

            if not self.sequence_file:
                self.steps = all_steps

            logger.info('Parsing AC Transport measurement.')
            logger.info(runs_list)
            self.data = split_ppms_data_act(data_df, runs_list)

        super().normalize(archive, logger)


m_package_ppms_act_default.__init_metainfo__()


m_package_ppms_mpms_default = SchemaPackage()


class CPFSPPMSMPMSMeasurementDefault(CPFSPPMSMPMSMeasurement, PlotSection, EntryData):
    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            # For automatic step discovery, some parameters are needed:
            if not self.temperature_tolerance:
                self.temperature_tolerance = 0.2
            if not self.field_tolerance:
                self.field_tolerance = 5.0

            with archive.m_context.raw_file(self.data_file, 'r') as file:
                data = file.read()

            header_match = re.search(r'\[Header\](.*?)\[Data\]', data, re.DOTALL)
            header_section = header_match.group(1).strip()
            header_lines = header_section.split('\n')

            while self.samples:
                self.m_remove_sub_section(CPFSPPMSMPMSMeasurementDefault.samples, 0)

            sample_headers = [
                line
                for line in header_lines
                if line.startswith('INFO') and 'SAMPLE_' in line
            ]
            if sample_headers:
                sample = CPFSSample()
                for line in sample_headers:
                    parts = re.split(r',\s*', line)
                    key = parts[-1].lower().replace('sample_', '')
                    if key == 'comment':
                        setattr(sample, key, ', '.join(parts[1:-1]))
                        # ids="_".join(parts[1].split("_")[1:3])
                        # ids = parts[1].split('_')[1]
                        # logger.info(ids)
                        # search_result = search(
                        #     owner='user',
                        #     query={
                        #         'results.eln.sections:any': ['CPFSCrystal'],
                        #         'results.eln.names:any': [ids + r'*'],
                        #     },
                        #     user_id=archive.metadata.main_author.user_id,
                        # )
                        # if len(search_result.data) > 0:
                        #     sample.sample_id = f"../uploads/{search_result.data[0]['upload_id']}/archive/{search_result.data[0]['entry_id']}#data"
                        #     sample.name = search_result.data[0]['search_quantities'][0][
                        #         'str_value'
                        #     ]
                        # else:
                        #     logger.warning(
                        #         "The sample given in the header could not be found and couldn't be referenced."
                        #     )
                    elif hasattr(sample, key):
                        setattr(sample, key, ', '.join(parts[1:-1]))
                self.m_add_sub_section(CPFSPPMSMPMSMeasurementDefault.samples, sample)

            for line in header_lines:
                if line.startswith('FILEOPENTIME'):
                    if hasattr(self, 'datetime'):
                        iso_date = get_fileopentime(line)
                        if iso_date == 'Not found.':
                            logger.error(
                                'FILEOPENTIME not understood. Check the format.'
                            )
                        else:
                            setattr(self, 'datetime', iso_date)
                if line.startswith('BYAPP'):
                    if hasattr(self, 'software'):
                        setattr(self, 'software', line.replace('BYAPP,', '').strip())
                if line.startswith('TEMPERATURETOLERANCE'):
                    if hasattr(self, 'temperature_tolerance'):
                        setattr(
                            self,
                            'temperature_tolerance',
                            float(line.replace('TEMPERATURETOLERANCE,', '').strip()),
                        )
                if line.startswith('FIELDTOLERANCE'):
                    if hasattr(self, 'field_tolerance'):
                        setattr(
                            self,
                            'field_tolerance',
                            float(line.replace('FIELDTOLERANCE,', '').strip()),
                        )

            data_section = header_match.string[header_match.end() :]
            data_section = data_section.replace(',Field', ',Magnetic Field')
            data_buffer = StringIO(data_section)
            data_df = pd.read_csv(
                data_buffer,
                header=0,
                skipinitialspace=True,
                sep=',',
                engine='python',
            )

            all_steps, runs_list = get_ppms_steps_from_data(
                data_df, self.temperature_tolerance, self.field_tolerance
            )

            if not self.sequence_file:
                self.steps = all_steps

            logger.info('Parsing MPMS measurement.')
            logger.info(runs_list)
            self.data = split_ppms_data_mpms(data_df, runs_list)

        super().normalize(archive, logger)


m_package_ppms_mpms_default.__init_metainfo__()


m_package_ppms_acms_default = SchemaPackage()


class CPFSPPMSACMSMeasurementDefault(CPFSPPMSACMSMeasurement, PlotSection, EntryData):
    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            # For automatic step discovery, some parameters are needed:
            if not self.temperature_tolerance:
                self.temperature_tolerance = 0.2
            if not self.field_tolerance:
                self.field_tolerance = 5.0

            with archive.m_context.raw_file(self.data_file, 'r') as file:
                data = file.read()

            header_match = re.search(r'\[Header\](.*?)\[Data\]', data, re.DOTALL)
            header_section = header_match.group(1).strip()
            header_lines = header_section.split('\n')

            while self.samples:
                self.m_remove_sub_section(CPFSPPMSACMSMeasurementDefault.samples, 0)

            sample_headers = [
                line
                for line in header_lines
                if line.startswith('INFO') and 'SAMPLE_' in line
            ]
            if sample_headers:
                sample = CPFSSample()
                for line in sample_headers:
                    parts = re.split(r',\s*', line)
                    key = parts[-1].lower().replace('sample_', '')
                    if key == 'comment':
                        setattr(sample, key, ', '.join(parts[1:-1]))
                        # ids="_".join(parts[1].split("_")[1:3])
                        # ids = parts[1].split('_')[1]
                        # logger.info(ids)
                        # search_result = search(
                        #     owner='user',
                        #     query={
                        #         'results.eln.sections:any': ['CPFSCrystal'],
                        #         'results.eln.names:any': [ids + r'*'],
                        #     },
                        #     user_id=archive.metadata.main_author.user_id,
                        # )
                        # if len(search_result.data) > 0:
                        #     sample.sample_id = f"../uploads/{search_result.data[0]['upload_id']}/archive/{search_result.data[0]['entry_id']}#data"
                        #     sample.name = search_result.data[0]['search_quantities'][0][
                        #         'str_value'
                        #     ]
                        # else:
                        #     logger.warning(
                        #         "The sample given in the header could not be found and couldn't be referenced."
                        #     )
                    elif hasattr(sample, key):
                        setattr(sample, key, ', '.join(parts[1:-1]))
                self.m_add_sub_section(CPFSPPMSACMSMeasurementDefault.samples, sample)

            for line in header_lines:
                if line.startswith('FILEOPENTIME'):
                    if hasattr(self, 'datetime'):
                        iso_date = get_fileopentime(line)
                        if iso_date == 'Not found.':
                            logger.error(
                                'FILEOPENTIME not understood. Check the format.'
                            )
                        else:
                            setattr(self, 'datetime', iso_date)
                if line.startswith('BYAPP'):
                    if hasattr(self, 'software'):
                        setattr(self, 'software', line.replace('BYAPP,', '').strip())
                if line.startswith('TEMPERATURETOLERANCE'):
                    if hasattr(self, 'temperature_tolerance'):
                        setattr(
                            self,
                            'temperature_tolerance',
                            float(line.replace('TEMPERATURETOLERANCE,', '').strip()),
                        )
                if line.startswith('FIELDTOLERANCE'):
                    if hasattr(self, 'field_tolerance'):
                        setattr(
                            self,
                            'field_tolerance',
                            float(line.replace('FIELDTOLERANCE,', '').strip()),
                        )

            data_section = header_match.string[header_match.end() :]
            data_section = data_section.replace(',Field', ',Magnetic Field')
            data_buffer = StringIO(data_section)
            data_df = pd.read_csv(
                data_buffer,
                header=0,
                skipinitialspace=True,
                sep=',',
                engine='python',
            )

            all_steps, runs_list = get_ppms_steps_from_data(
                data_df, self.temperature_tolerance, self.field_tolerance
            )

            if not self.sequence_file:
                self.steps = all_steps

            logger.info('Parsing ACMS measurement.')
            logger.info(runs_list)
            self.data = split_ppms_data_acms(data_df, runs_list)

        super().normalize(archive, logger)


m_package_ppms_acms_default.__init_metainfo__()
