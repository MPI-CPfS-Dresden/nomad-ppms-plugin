from nomad.config.models.plugins import (
    ParserEntryPoint,
    SchemaPackageEntryPoint,
)
from pydantic import Field


class CPFSPPMSETOEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_eto_default,
        )

        return m_package_ppms_eto_default


class CPFSPPMSETOEntryPointLabview(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_eto_labview,
        )

        return m_package_ppms_eto_labview


class CPFSPPMSACTEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_act_default,
        )

        return m_package_ppms_act_default


class CPFSPPMSMPMSEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_mpms_default,
        )

        return m_package_ppms_mpms_default


class CPFSPPMSACMSEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_acms_default,
        )

        return m_package_ppms_acms_default


class CPFSPPMSResistivityEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.schema_package import (
            m_package_ppms_resistivity_default,
        )

        return m_package_ppms_resistivity_default


eto_schema_default = CPFSPPMSETOEntryPointDefault(
    name='CPFSPPMSETOEntryPoint',
    description='New schema package entry point configuration.',
)

eto_schema_labview = CPFSPPMSETOEntryPointLabview(
    name='CPFSPPMSETOEntryPoint',
    description='New schema package entry point configuration.',
)

act_schema_default = CPFSPPMSACTEntryPointDefault(
    name='CPFSPPMSACTEntryPoint',
    description='New schema package entry point configuration.',
)

mpms_schema_default = CPFSPPMSMPMSEntryPointDefault(
    name='CPFSPPMSMPMSEntryPoint',
    description='New schema package entry point configuration.',
)

acms_schema_default = CPFSPPMSACMSEntryPointDefault(
    name='CPFSPPMSACMSEntryPoint',
    description='New schema package entry point configuration.',
)
resistivity_schema_default = CPFSPPMSResistivityEntryPointDefault(
    name='CPFSPPMSACMSEntryPoint',
    description='New schema package entry point configuration.',
)


class DataParserEntryPointETODefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSETOParserDefault

        return CPFSPPMSETOParserDefault(**self.dict())


eto_parser_default = DataParserEntryPointETODefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain|application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP, Electrical Transport Option',
)


class DataParserEntryPointETOLabview(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSETOParserLabview

        return CPFSPPMSETOParserLabview(**self.dict())


eto_parser_labview = DataParserEntryPointETOLabview(
    name='DataParser',
    description='New parser entry point configuration.',
    # mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain',
    mainfile_contents_re='; LABVIEW measurement file V2.6',
)


class DataParserEntryPointACTDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSACTParserDefault

        return CPFSPPMSACTParserDefault(**self.dict())


act_parser_default = DataParserEntryPointACTDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain|application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*ACTRANSPORT',
)


class DataParserEntryPointMPMSDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSMPMSParserDefault

        return CPFSPPMSMPMSParserDefault(**self.dict())


mpms_parser_default = DataParserEntryPointMPMSDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain|application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*MPMS3',
)


class DataParserEntryPointACMSDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSACMSParserDefault

        return CPFSPPMSACMSParserDefault(**self.dict())


acms_parser_default = DataParserEntryPointACMSDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain|application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*ACMS',
)


class DataParserEntryPointResistivityDefault(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_ppms_plugin.parser import CPFSPPMSResisitivityParserDefault

        return CPFSPPMSResisitivityParserDefault(**self.dict())


resistivity_parser_default = DataParserEntryPointResistivityDefault(
    name='DataParser',
    description='New parser entry point configuration.',
    mainfile_name_re='^.+\.dat$',
    mainfile_mime_re='text/plain|application/x-wine-extension-ini',
    mainfile_contents_re='BYAPP,\s*Resistivity',
)
