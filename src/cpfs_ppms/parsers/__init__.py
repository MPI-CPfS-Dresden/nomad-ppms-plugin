from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class DataParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.parsers.parser import CPFSPPMSParser

        return CPFSPPMSParser(**self.dict())


parser_entry_point_data = DataParserEntryPoint(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='application/x-wine-extension-ini',
)


class SqcParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.parsers.parser import CPFSPPMSSequenceParser

        return CPFSPPMSSequenceParser(**self.dict())


parser_entry_point_sqc = SqcParserEntryPoint(
    name='SequenceParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.seq$',
    mainfile_mime_re='text/plain',
)
