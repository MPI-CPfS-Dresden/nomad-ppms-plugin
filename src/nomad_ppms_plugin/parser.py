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

from typing import (
    TYPE_CHECKING,
)

from nomad.datamodel import EntryArchive
from nomad_measurements.ppms.parser import (
    PPMSACMSParser,
    PPMSACTParser,
    PPMSETOParser,
    PPMSMPMSParser,
    PPMSResistivityParser,
)

from nomad_ppms_plugin.schema_package import (
    CPFSPPMSACMSMeasurementDefault,
    CPFSPPMSACTMeasurementDefault,
    CPFSPPMSETOMeasurementDefault,
    CPFSPPMSETOMeasurementLabview,
    CPFSPPMSMPMSMeasurementDefault,
    CPFSPPMSResistivityMeasurementDefault,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )

from nomad.config import config
from nomad.datamodel import EntryArchive

configuration = config.get_plugin_entry_point('nomad_ppms_plugin:eto_parser_default')
configuration = config.get_plugin_entry_point('nomad_ppms_plugin:eto_parser_labview')
configuration = config.get_plugin_entry_point('nomad_ppms_plugin:act_parser_default')
configuration = config.get_plugin_entry_point('nomad_ppms_plugin:mpms_parser_default')
configuration = config.get_plugin_entry_point('nomad_ppms_plugin:acms_parser_default')
configuration = config.get_plugin_entry_point(
    'nomad_ppms_plugin:resistivity_parser_default'
)


class CPFSPPMSETOParserDefault(PPMSETOParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSETOMeasurementDefault

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)


class CPFSPPMSETOParserLabview(PPMSETOParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSETOMeasurementLabview

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)


class CPFSPPMSACTParserDefault(PPMSACTParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSACTMeasurementDefault

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)


class CPFSPPMSMPMSParserDefault(PPMSMPMSParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSMPMSMeasurementDefault

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)


class CPFSPPMSACMSParserDefault(PPMSACMSParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSACMSMeasurementDefault

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)


class CPFSPPMSResisitivityParserDefault(PPMSResistivityParser):
    def set_entrydata_definition(self):
        self.entrydata_definition = CPFSPPMSResistivityMeasurementDefault

    def parse(self, mainfile: str, archive: EntryArchive, logger) -> None:
        super().parse(mainfile, archive, logger)
