from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class CPFSPPMSETOEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_eto_default

        return m_package_ppms_eto_default


class CPFSPPMSETOEntryPointLabview(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_eto_labview

        return m_package_ppms_eto_labview


class CPFSPPMSACTEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_act_default

        return m_package_ppms_act_default


class CPFSPPMSMPMSEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_mpms_default

        return m_package_ppms_mpms_default


class CPFSPPMSACMSEntryPointDefault(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_acms_default

        return m_package_ppms_acms_default


schema_entry_point_eto_default = CPFSPPMSETOEntryPointDefault(
    name='CPFSPPMSETOEntryPoint',
    description='New schema package entry point configuration.',
)

schema_entry_point_eto_labview = CPFSPPMSETOEntryPointLabview(
    name='CPFSPPMSETOEntryPoint',
    description='New schema package entry point configuration.',
)

schema_entry_point_act_default = CPFSPPMSACTEntryPointDefault(
    name='CPFSPPMSACTEntryPoint',
    description='New schema package entry point configuration.',
)

schema_entry_point_mpms_default = CPFSPPMSMPMSEntryPointDefault(
    name='CPFSPPMSMPMSEntryPoint',
    description='New schema package entry point configuration.',
)

schema_entry_point_acms_default = CPFSPPMSACMSEntryPointDefault(
    name='CPFSPPMSACMSEntryPoint',
    description='New schema package entry point configuration.',
)
