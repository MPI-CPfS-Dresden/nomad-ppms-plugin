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

from time import perf_counter, sleep
from typing import (
    TYPE_CHECKING,
)

from nomad.datamodel import EntryArchive
from nomad.datamodel.data import (
    EntryData,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
)
from nomad.metainfo import Quantity
from nomad.parsing import MatchingParser
from nomad.search import search
from nomad_material_processing.utils import create_archive

from nomad_ppms_plugin.schema_packages.schema_package import (
    CPFSPPMSACMSMeasurementDefault,
    CPFSPPMSACTMeasurementDefault,
    CPFSPPMSETOMeasurementDefault,
    CPFSPPMSETOMeasurementLabview,
    CPFSPPMSMPMSMeasurementDefault,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )

from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.basesections import (
    BaseSection,
)

configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_data_eto_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_data_eto_labview'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_data_act_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_data_mpms_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_data_acms_default'
)
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin.parsers:parser_entry_point_sqc'
)


class CPFSPPMSFile(EntryData):
    measurement = Quantity(
        type=CPFSPPMSETOMeasurementDefault,
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
        ),
    )


class CPFSPPMSETOParserDefault(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        timethreshold = 15
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        entry = CPFSPPMSETOMeasurementDefault()
        entry.data_file = data_file_with_path
        file_name = f'{data_file[:-4]}.archive.json'
        # entry.normalize(archive, logger)
        tic = perf_counter()
        while True:
            search_result = search(
                owner='user',
                query={
                    'results.eln.sections:any': ['CPFSPPMSSequenceFile'],
                    'upload_id:any': [archive.m_context.upload_id],
                },
                user_id=archive.metadata.main_author.user_id,
            )
            if len(search_result.data) > 0:
                for sequence in search_result.data:
                    entry.sequence_file = sequence['search_quantities'][0]['str_value']
                    logger.info(sequence['search_quantities'][0]['str_value'])
                    break
            sleep(0.1)
            toc = perf_counter()
            if toc - tic > timethreshold:
                logger.warning(
                    "The Sequence File entry/ies in the current upload \
                               were not found and couldn't be referenced."
                )
                break
        archive.data = CPFSPPMSFile(
            measurement=create_archive(entry, archive, file_name)
        )
        archive.metadata.entry_name = data_file + ' measurement file'


class CPFSPPMSETOParserLabview(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        entry = CPFSPPMSETOMeasurementLabview()
        entry.data_file = data_file_with_path
        file_name = f'{data_file[:-4]}.archive.json'
        archive.data = CPFSPPMSFile(
            measurement=create_archive(entry, archive, file_name)
        )
        archive.metadata.entry_name = data_file + ' measurement file'


class CPFSPPMSACTParserDefault(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        timethreshold = 15
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        entry = CPFSPPMSACTMeasurementDefault()
        entry.data_file = data_file_with_path
        file_name = f'{data_file[:-4]}.archive.json'
        # entry.normalize(archive, logger)
        tic = perf_counter()
        while True:
            search_result = search(
                owner='user',
                query={
                    'results.eln.sections:any': ['CPFSPPMSSequenceFile'],
                    'upload_id:any': [archive.m_context.upload_id],
                },
                user_id=archive.metadata.main_author.user_id,
            )
            if len(search_result.data) > 0:
                for sequence in search_result.data:
                    entry.sequence_file = sequence['search_quantities'][0]['str_value']
                    logger.info(sequence['search_quantities'][0]['str_value'])
                    break
            sleep(0.1)
            toc = perf_counter()
            if toc - tic > timethreshold:
                logger.warning(
                    "The Sequence File entry/ies in the current upload \
                               were not found and couldn't be referenced."
                )
                break
        archive.data = CPFSPPMSFile(
            measurement=create_archive(entry, archive, file_name)
        )
        archive.metadata.entry_name = data_file + ' measurement file'


class CPFSPPMSMPMSParserDefault(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        timethreshold = 15
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        entry = CPFSPPMSMPMSMeasurementDefault()
        entry.data_file = data_file_with_path
        file_name = f'{data_file[:-4]}.archive.json'
        # entry.normalize(archive, logger)
        tic = perf_counter()
        while True:
            search_result = search(
                owner='user',
                query={
                    'results.eln.sections:any': ['CPFSPPMSSequenceFile'],
                    'upload_id:any': [archive.m_context.upload_id],
                },
                user_id=archive.metadata.main_author.user_id,
            )
            if len(search_result.data) > 0:
                for sequence in search_result.data:
                    entry.sequence_file = sequence['search_quantities'][0]['str_value']
                    logger.info(sequence['search_quantities'][0]['str_value'])
                    break
            sleep(0.1)
            toc = perf_counter()
            if toc - tic > timethreshold:
                logger.warning(
                    "The Sequence File entry/ies in the current upload \
                               were not found and couldn't be referenced."
                )
                break
        archive.data = CPFSPPMSFile(
            measurement=create_archive(entry, archive, file_name)
        )
        archive.metadata.entry_name = data_file + ' measurement file'


class CPFSPPMSACMSParserDefault(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        timethreshold = 15
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        entry = CPFSPPMSACMSMeasurementDefault()
        entry.data_file = data_file_with_path
        file_name = f'{data_file[:-4]}.archive.json'
        # entry.normalize(archive, logger)
        tic = perf_counter()
        while True:
            search_result = search(
                owner='user',
                query={
                    'results.eln.sections:any': ['CPFSPPMSSequenceFile'],
                    'upload_id:any': [archive.m_context.upload_id],
                },
                user_id=archive.metadata.main_author.user_id,
            )
            if len(search_result.data) > 0:
                for sequence in search_result.data:
                    entry.sequence_file = sequence['search_quantities'][0]['str_value']
                    logger.info(sequence['search_quantities'][0]['str_value'])
                    break
            sleep(0.1)
            toc = perf_counter()
            if toc - tic > timethreshold:
                logger.warning(
                    "The Sequence File entry/ies in the current upload \
                               were not found and couldn't be referenced."
                )
                break
        archive.data = CPFSPPMSFile(
            measurement=create_archive(entry, archive, file_name)
        )
        archive.metadata.entry_name = data_file + ' measurement file'


class CPFSPPMSSequenceFile(BaseSection, EntryData):
    file_path = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
    )


class CPFSPPMSSequenceParser(MatchingParser):
    # def __init__(self):
    #    super().__init__(
    #        name='NOMAD PPMS schema and parser plugin for the CPFS',
    #        code_name='nomad_ppms_plugin_sequence',
    #        code_homepage='https://github.com/FAIRmat-NFDI/AreaA-data_modeling_and_schemas',
    #        supported_compressions=['gz', 'bz2', 'xz'],
    #    )

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        data_file = mainfile.split('/')[-1]
        data_file_with_path = mainfile.split('raw/')[-1]
        # file_name = f'{data_file[:-4]}.archive.json'
        # entry.normalize(archive, logger)
        archive.data = CPFSPPMSSequenceFile(file_path=data_file_with_path)
        archive.metadata.entry_name = data_file + ' sequence file'
