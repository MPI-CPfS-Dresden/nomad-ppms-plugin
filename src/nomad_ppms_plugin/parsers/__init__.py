from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class DataParserEntryPointETODefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSETOParserDefault

        return CPFSPPMSETOParserDefault(**self.dict())


parser_entry_point_data_eto_default = DataParserEntryPointETODefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP, Electrical Transport Option',
)


class DataParserEntryPointETOLabview(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSETOParserLabview

        return CPFSPPMSETOParserLabview(**self.dict())


parser_entry_point_data_eto_labview = DataParserEntryPointETOLabview(
    name='DataParser',
    description='New parser entry point configuration.',
    # mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain',
    mainfile_contents_re='; LABVIEW measurement file V2.6',
)


class DataParserEntryPointACTDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSACTParserDefault

        return CPFSPPMSACTParserDefault(**self.dict())


parser_entry_point_data_act_default = DataParserEntryPointACTDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*ACTRANSPORT',
)


class DataParserEntryPointMPMSDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSMPMSParserDefault

        return CPFSPPMSMPMSParserDefault(**self.dict())


parser_entry_point_data_mpms_default = DataParserEntryPointMPMSDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*MPMS3',
)


class DataParserEntryPointACMSDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSACMSParserDefault

        return CPFSPPMSACMSParserDefault(**self.dict())


parser_entry_point_data_acms_default = DataParserEntryPointACMSDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*ACMS',
)


class SqcParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parsers.parser import CPFSPPMSSequenceParser

        return CPFSPPMSSequenceParser(**self.dict())


parser_entry_point_sqc = SqcParserEntryPoint(
    name='SequenceParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.seq$',
    mainfile_mime_re='text/plain',
)
