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

import numpy as np
import pandas as pd
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    SectionProperties,
)
from nomad.datamodel.data import (
    EntryData,
)
from nomad.datamodel.metainfo.basesections import Measurement
from nomad.datamodel.metainfo.plot import PlotlyFigure, PlotSection
from nomad.metainfo import (
    MEnum,
    Quantity,
    Section,
    SubSection,
)
from structlog.stdlib import (
    BoundLogger,
)
from nomad.search import search

from nomad.units import ureg

from cpfs_ppms.ppmsdatastruct import (
    PPMSData,
)
from cpfs_ppms.ppmsfunctions import (
    find_ppms_steps_from_sequence,
    get_fileopentime,
    get_ppms_steps_from_data,
    split_ppms_data_act,
    split_ppms_data_eto,
)
from cpfs_ppms.ppmssteps import (
    PPMSMeasurementStep,
)
from cpfs_ppms.cpfsppmsdatastruct import (
    CPFSETOAnalyzedData,
    CPFSETOSymmetrizedData,
    CPFSSample,
    CPFSPPMSETOMeasurement,
    CPFSPPMSMeasurement,
)

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.config import config
from nomad.metainfo import SchemaPackage

configuration = config.get_plugin_entry_point(
    'cpfs_ppms.schema_packages:schema_entry_point_eto_default'
)

m_package_ppms_eto_default = SchemaPackage()


class CPFSPPMSETOMeasurementDefault(CPFSPPMSETOMeasurement, PlotSection, EntryData):


    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915

        if archive.data.data_file:
            logger.info('Parsing PPMS measurement file.')
            #For automatic step discovery, some parameters are needed:
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
                                sample.sample_id = f"../uploads/{search_result.data[0]['upload_id']}/archive/{search_result.data[0]['entry_id']}#data"
                                sample.name = search_result.data[0]['search_quantities'][0][
                                    'str_value'
                                ]
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
                    self.m_add_sub_section(CPFSPPMSETOMeasurementDefault.samples, sample)

            for line in header_lines:
                if line.startswith('FILEOPENTIME'):
                    if hasattr(self, 'datetime'):
                        iso_date = get_fileopentime(line)
                        if iso_date=="Not found.":
                            logger.error("FILEOPENTIME not understood. Check the format.")
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
                data_df,self.temperature_tolerance,self.field_tolerance
                )

            if not self.sequence_file:
                self.steps=all_steps

            if self.software.startswith('ACTRANSPORT'):
                logger.info('Parsing AC Transport measurement.')
                logger.info(runs_list)
                self.data=split_ppms_data_act(data_df,runs_list)

            if self.software.startswith('Electrical Transport Option'):
                logger.info('Parsing ETO measurement.')
                logger.info(runs_list)
                self.data=split_ppms_data_eto(data_df,runs_list)

        super().normalize(archive, logger)


m_package_ppms_eto_default.__init_metainfo__()
