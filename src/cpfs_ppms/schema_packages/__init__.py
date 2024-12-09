from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class CPFSPPMSETOEntryPoint(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from cpfs_ppms.schema_packages.schema_package import m_package_ppms_eto_default

        return m_package_ppms_eto_default


schema_entry_point_eto_default = CPFSPPMSETOEntryPoint(
    name='CPFSPPMSETOEntryPoint',
    description='New schema package entry point configuration.',
)
