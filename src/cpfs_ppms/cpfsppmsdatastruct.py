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
import pandas as pd
from nomad.datamodel.data import (
    EntryData,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    SectionProperties,
)
from nomad.datamodel.metainfo.basesections import Measurement
from nomad.datamodel.metainfo.eln import (
    CompositeSystem,
    SampleID,
)
from nomad.datamodel.metainfo.plot import PlotlyFigure, PlotSection
from nomad.metainfo import (
    Datetime,
    MEnum,
    Quantity,
    Section,
    SubSection,
)
from nomad.units import ureg
from nomad_ppms_plugin.ppmsdatastruct import PPMSData
from nomad_ppms_plugin.ppmsfunctions import (
    find_ppms_steps_from_sequence,
)
from nomad_ppms_plugin.ppmssteps import (
    PPMSMeasurementStep,
)
from structlog.stdlib import (
    BoundLogger,
)


class CPFSCrystal(EntryData):
    sample_id = Quantity(
        type=str,
        description="""
        Sample ID given by the grower.
        """,
        a_eln={
            'component': 'StringEditQuantity',
        },
    )

    internal_sample_id = SubSection(
        section_def=SampleID,
    )
    achieved_composition = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    final_crystal_length = Quantity(
        type=float,
        unit='meter',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity',
            defaultDisplayUnit='millimeter',
        ),
    )
    single_poly = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    crystal_shape = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    crystal_orientation = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    safety_reactivity = Quantity(
        type=str,
        a_eln=ELNAnnotation(
            component='StringEditQuantity',
        ),
    )
    datetime = Quantity(
        type=Datetime,
        description='The date and time associated with this section.',
    )
    lab_id = Quantity(
        type=str,
        description="""An ID string that is unique at least for the lab that produced this
            data.""",  # noqa: E501
    )
    description = Quantity(
        type=str,
        description='Any information that cannot be captured in the other fields.',
        a_eln=dict(component='StringEditQuantity', label='Remarks'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:
        """
        The normalizer for the `CPFSCrystal` class.

        Args:
            archive (EntryArchive): The archive containing the section that is being
            normalized.
            logger (BoundLogger): A structlog logger.
        """
        super().normalize(archive, logger)


class CPFSSample(CompositeSystem):
    name = Quantity(type=str, description='FILL')
    type = Quantity(type=str, description='FILL')
    material = Quantity(type=str, description='FILL')
    comment = Quantity(type=str, description='FILL')
    lead_separation = Quantity(type=str, description='FILL')
    cross_section = Quantity(type=str, description='FILL')
    sample_id = Quantity(
        type=CPFSCrystal,
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
        ),
    )
    width = Quantity(
        type=np.float64,
        unit='meter',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='millimeter'
        ),
    )
    length = Quantity(
        type=np.float64,
        unit='meter',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='millimeter'
        ),
    )
    depth = Quantity(
        type=np.float64,
        unit='meter',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='millimeter'
        ),
    )


class CPFSETOSymmetrizedData(PPMSData):
    field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    rho_xx_up = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_xx_down = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    mr_up = Quantity(
        type=np.dtype(np.float64), shape=['*'], unit='dimensionless', description='FILL'
    )
    mr_down = Quantity(
        type=np.dtype(np.float64), shape=['*'], unit='dimensionless', description='FILL'
    )
    rho_xy_up = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_xy_down = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )


class CPFSETOAnalyzedData(PPMSData):
    field = Quantity(
        type=np.dtype(np.float64), unit='gauss', shape=['*'], description='FILL'
    )
    rho_xx_up = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_xx_down = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    sigma_xx_up = Quantity(
        type=np.dtype(np.float64), unit='siemens/m', shape=['*'], description='FILL'
    )
    sigma_xx_down = Quantity(
        type=np.dtype(np.float64), unit='siemens/m', shape=['*'], description='FILL'
    )
    rho_ohe_up = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_ohe_down = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_ahe_up = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    rho_ahe_down = Quantity(
        type=np.dtype(np.float64), unit='ohm*meter', shape=['*'], description='FILL'
    )
    sigma_ahe_up = Quantity(
        type=np.dtype(np.float64), unit='siemens/m', shape=['*'], description='FILL'
    )
    sigma_ahe_down = Quantity(
        type=np.dtype(np.float64), unit='siemens/m', shape=['*'], description='FILL'
    )
    carrier_concentration = Quantity(
        type=np.dtype(np.float64), unit='1/meter**3', description='FILL'
    )
    carrier_mobility = Quantity(
        type=np.dtype(np.float64), unit='meter**2/volt/second', description='FILL'
    )


class CPFSPPMSMeasurement(Measurement):
    m_def = Section(
        a_eln=ELNAnnotation(
            properties=SectionProperties(
                order=[
                    'name',
                    'datetime',
                    'data_file',
                    'sequence_file',
                    'description',
                    'software',
                    'startupaxis',
                ],
            ),
            lane_width='600px',
        ),
    )

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
    )
    file_open_time = Quantity(type=str, description='FILL')
    software = Quantity(type=str, description='FILL')

    steps = SubSection(
        section_def=PPMSMeasurementStep,
        repeats=True,
    )

    data = SubSection(section_def=PPMSData, repeats=True)

    sequence_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        pass
        # super().normalize(archive, logger)


class CPFSPPMSETOMeasurement(CPFSPPMSMeasurement, PlotSection, EntryData):
    temperature_tolerance = Quantity(
        type=float,
        unit='kelvin',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='kelvin'
        ),
    )

    field_tolerance = Quantity(
        type=float,
        unit='gauss',
        a_eln=ELNAnnotation(component='NumberEditQuantity', defaultDisplayUnit='gauss'),
    )

    symmetrized_data = SubSection(
        section_def=CPFSETOSymmetrizedData,
        repeats=True,
    )

    analyzed_data = SubSection(
        section_def=CPFSETOAnalyzedData,
        repeats=True,
    )

    channel_measurement_type = Quantity(
        type=MEnum(
            'TMR',
            'Hall',
            'undefined',
        ),
        shape=['*'],
        a_eln=ELNAnnotation(
            component='EnumEditQuantity',
        ),
    )

    testfield2 = Quantity(
        type=float,
    )

    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        super().normalize(archive, logger)

        ### Start of the PPMSMeasurement normalizer

        logger.info('AAA')

        self.testfield = 1.0

        if archive.data.sequence_file:
            logger.info('Parsing PPMS sequence file.')
            with archive.m_context.raw_file(self.sequence_file, 'r') as file:
                sequence = file.readlines()
                self.steps = find_ppms_steps_from_sequence(sequence)

        # Now create the according plots
        import plotly.express as px
        from plotly.subplots import make_subplots

        for data in self.data:
            if data.measurement_type == 'field':
                resistivity_ch1 = px.scatter(
                    x=data.magnetic_field, y=data.channels[0].resistance
                )
                resistivity_ch2 = px.scatter(
                    x=data.magnetic_field, y=data.channels[1].resistance
                )
            if data.measurement_type == 'temperature':
                resistivity_ch1 = px.scatter(
                    x=data.temperature, y=data.channels[0].resistance
                )
                resistivity_ch2 = px.scatter(
                    x=data.temperature, y=data.channels[1].resistance
                )
            figure1 = make_subplots(rows=2, cols=1, shared_xaxes=True)
            figure1.add_trace(resistivity_ch1.data[0], row=1, col=1)
            figure1.add_trace(resistivity_ch2.data[0], row=2, col=1)
            figure1.update_layout(height=400, width=716, title_text=data.name)
            self.figures.append(
                PlotlyFigure(label=data.name, figure=figure1.to_plotly_json())
            )

        # find measurement modes, for now coming from sample.comment
        if self.samples[0].comment:
            modelist = []
            for channel in ['Ch1_', 'Ch2_']:
                if channel + 'TMR' in self.samples[0].comment:
                    modelist.append('TMR')
                elif channel + 'Hall' in self.samples[0].comment:
                    modelist.append('Hall')
                else:
                    modelist.append('undefined')
            self.channel_measurement_type = modelist

        if (
            'Hall' in self.channel_measurement_type
            and 'TMR' in self.channel_measurement_type
        ):
            # find biggest fitlength
            maxfield = 90000 * ureg('gauss')
            fitlength = 0
            for mdata in self.data:
                fitlength = max(len(mdata.magnetic_field) / 4, fitlength)
            fitlength -= fitlength % -100
            fitlength += 1
            fitfield = np.linspace(-maxfield, maxfield, int(fitlength))
            # Try to symmetrize data for each measurement
            data_symmetrized = []
            for mdata in self.data:
                # For now only for field sweeps
                if not mdata.name.startswith('Field sweep'):
                    continue
                sym_data = CPFSETOSymmetrizedData()
                sym_data.name = mdata.name
                sym_data.title = mdata.name
                sym_data.field = fitfield
                for channel in [0, 1]:
                    field = mdata.magnetic_field[
                        np.invert(pd.isnull(mdata.channels[channel].resistance))
                    ]
                    res = mdata.channels[channel].resistance[
                        np.invert(pd.isnull(mdata.channels[channel].resistance))
                    ]
                    # Check if field sweeps down and up:
                    downsweep = []
                    upsweep = []
                    for i in range(len(field)):
                        if abs(field[i] - maxfield) < self.field_tolerance:
                            if len(downsweep) == 0:
                                downsweep.append(i)  # downsweep started
                            if len(upsweep) == 1:
                                upsweep.append(i)  # upsweep finished
                        if abs(field[i] + maxfield) < self.field_tolerance:
                            if len(downsweep) == 1:
                                downsweep.append(i)  # downsweep finished
                            if len(upsweep) == 0:
                                upsweep.append(i)  # upsweep started
                    if len(upsweep) != 2 and len(downsweep) != 2:  # noqa: PLR2004
                        logger.warning(
                            'Measurement '
                            + mdata.name
                            + ' did not contain up- and downsweep in field.'
                        )
                        continue
                    elif len(downsweep) != 2:  # noqa: PLR2004
                        logger.warning(
                            'Measurement '
                            + mdata.name
                            + ' did not contain downsweep in field.\
                                  Using upsweep for both ways.'
                        )
                        upfit = np.interp(
                            fitfield,
                            field[upsweep[0] : upsweep[1]],
                            res[upsweep[0] : upsweep[1]],
                        )
                        downfit = upfit
                    elif len(upsweep) != 2:  # noqa: PLR2004
                        logger.warning(
                            'Measurement '
                            + mdata.name
                            + ' did not contain upsweep in field.\
                                  Using downsweep for both ways.'
                        )
                        downfit = np.interp(
                            fitfield,
                            np.flip(field[downsweep[0] : downsweep[1]]),
                            np.flip(res[downsweep[0] : downsweep[1]]),
                        )
                        upfit = downfit
                    else:
                        downfit = np.interp(
                            fitfield,
                            np.flip(field[downsweep[0] : downsweep[1]]),
                            np.flip(res[downsweep[0] : downsweep[1]]),
                        )
                        upfit = np.interp(
                            fitfield,
                            field[upsweep[0] : upsweep[1]],
                            res[upsweep[0] : upsweep[1]],
                        )
                    if self.channel_measurement_type[channel] == 'Hall':
                        intermediate = (upfit + np.flip(downfit)) / 2.0
                        sym_data.rho_xy_down = (
                            downfit - np.flip(intermediate)
                        ) * self.samples[channel].depth
                        sym_data.rho_xy_up = (upfit - intermediate) * self.samples[
                            channel
                        ].depth
                    if self.channel_measurement_type[channel] == 'TMR':
                        intermediate = (np.flip(downfit) - upfit) / 2.0
                        sym_data.rho_xx_down = (
                            (downfit + np.flip(intermediate))
                            * self.samples[channel].depth
                            * self.samples[channel].width
                            / self.samples[channel].length
                        )
                        sym_data.rho_xx_up = (
                            (upfit - intermediate)
                            * self.samples[channel].depth
                            * self.samples[channel].width
                            / self.samples[channel].length
                        )
                        sym_data.mr_down = (
                            sym_data.rho_xx_down
                            - sym_data.rho_xx_down[int(fitlength / 2)]
                        ) / sym_data.rho_xx_down[int(fitlength / 2)]
                        sym_data.mr_up = (
                            sym_data.rho_xx_up - sym_data.rho_xx_up[int(fitlength / 2)]
                        ) / sym_data.rho_xx_up[int(fitlength / 2)]
                data_symmetrized.append(sym_data)

                # create symmetrized output files
                filename = (
                    'symmetrized_data_'
                    + '_'.join(sym_data.name.split())
                    + '_'
                    + self.data_file.strip('.dat')
                )
                with archive.m_context.raw_file(filename, 'w') as outfile:
                    outfile.write(
                        '#Field (Oe)     rho_xx_up       rho_xx_down      mr_up  \
                                 mr_down         rho_xy_up       rho_xy_down      \n'
                    )
                    for i in range(int(fitlength)):
                        outfile.write(f'{fitfield[i].magnitude:16.8e}')
                        if sym_data.rho_xx_up is not None:
                            outfile.write(f'{sym_data.rho_xx_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if sym_data.rho_xx_down is not None:
                            outfile.write(f'{sym_data.rho_xx_down[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if sym_data.mr_up is not None:
                            outfile.write(f'{sym_data.mr_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if sym_data.mr_down is not None:
                            outfile.write(f'{sym_data.mr_down[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if sym_data.rho_xy_up is not None:
                            outfile.write(f'{sym_data.rho_xy_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if sym_data.rho_xy_down is not None:
                            outfile.write(f'{sym_data.rho_xy_down[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        outfile.write('\n')

            self.symmetrized_data = data_symmetrized

            # Analyze data: split in ordinary and anomalous, calculate conductivities
            filename = (
                'analyzed_data_carrier_mobility_and_concentration_'
                + '_'
                + self.data_file.strip('.dat')
            )
            carrierout = archive.m_context.raw_file(filename, 'w')
            carrierout.write(
                '#Temperature      carrier concentration    carrier mobility\n'
            )
            cutofffield = 50000 * ureg('gauss')
            data_analyzed = []
            for data in self.symmetrized_data:
                ana_data = CPFSETOAnalyzedData()
                ana_data.name = data.name
                ana_data.title = data.name
                ana_data.field = data.field
                ana_data.rho_xx_up = data.rho_xx_up
                ana_data.rho_xx_down = data.rho_xx_down
                ana_data.sigma_xx_up = 1.0 / data.rho_xx_up
                ana_data.sigma_xx_down = 1.0 / data.rho_xx_down

                fitstart = int(
                    len(data.field) * (180000 - cutofffield.magnitude) / 180000
                )
                rho_xy_up_fit = (
                    np.poly1d(
                        [
                            np.polyfit(
                                data.field[fitstart:].magnitude,
                                data.rho_xy_up[fitstart:].magnitude,
                                1,
                            )[0],
                            0,
                        ]
                    )
                    * data.rho_xy_up.units
                )
                rho_xy_down_fit = (
                    np.poly1d(
                        [
                            np.polyfit(
                                data.field[fitstart:].magnitude,
                                data.rho_xy_down[fitstart:].magnitude,
                                1,
                            )[0],
                            0,
                        ]
                    )
                    * data.rho_xy_down.units
                )
                ana_data.rho_ohe_up = rho_xy_up_fit(ana_data.field.magnitude)
                ana_data.rho_ohe_down = rho_xy_down_fit(ana_data.field.magnitude)
                ana_data.rho_ahe_up = data.rho_xy_up - ana_data.rho_ohe_up
                ana_data.rho_ahe_down = data.rho_xy_down - ana_data.rho_ohe_down
                ana_data.sigma_ahe_up = ana_data.rho_ahe_up / (
                    ana_data.rho_ahe_up**2 + data.rho_xx_up**2
                )
                ana_data.sigma_ahe_down = ana_data.rho_ahe_down / (
                    ana_data.rho_ahe_down**2 + data.rho_xx_down**2
                )

                ana_data.carrier_concentration = 1.0 / (
                    (
                        np.polyfit(
                            data.field[fitstart:].magnitude,
                            data.rho_xy_up[fitstart:].magnitude,
                            1,
                        )[0]
                        * data.rho_xy_up.units
                        / data.field.units
                    )
                    * (1.60217663 * 10**-19 * ureg('coulomb'))
                )
                ana_data.carrier_mobility = 1.0 / (
                    ana_data.carrier_concentration
                    * (1.60217663 * 10**-19 * ureg('coulomb'))
                    * ana_data.rho_xx_up[int(len(ana_data.rho_xx_up) / 2)]
                )

                data_analyzed.append(ana_data)

                # create analyzed output files
                carrierout.write(
                    f'{ana_data.name.split()[3]}       \
                          {ana_data.carrier_concentration.magnitude / 1000000.0}    \
                                {ana_data.carrier_mobility.magnitude * 10000.0}\n'
                )
                filename = (
                    'analyzed_data_'
                    + '_'.join(ana_data.name.split())
                    + '_'
                    + self.data_file.strip('.dat')
                )
                with archive.m_context.raw_file(filename, 'w') as outfile:
                    outfile.write(
                        '#Field (Oe)     rho_ohe_up      rho_ahe_up      rho_ahe_down \
                              sigma_ahe_up       sigma_ahe_down      \n'
                    )
                    for i in range(int(fitlength)):
                        outfile.write(f'{fitfield[i].magnitude:16.8e}')
                        if ana_data.rho_ohe_up is not None:
                            outfile.write(f'{ana_data.rho_ohe_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if ana_data.rho_ahe_up is not None:
                            outfile.write(f'{ana_data.rho_ahe_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if ana_data.rho_ahe_down is not None:
                            outfile.write(f'{ana_data.rho_ahe_down[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if ana_data.sigma_ahe_up is not None:
                            outfile.write(f'{ana_data.sigma_ahe_up[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        if ana_data.rho_ahe_down is not None:
                            outfile.write(f'{ana_data.rho_ahe_down[i].magnitude:16.8e}')
                        else:
                            outfile.write('NaN             ')
                        outfile.write('\n')

            self.analyzed_data = data_analyzed

            # Now create the according plots
            import plotly.express as px
            import plotly.graph_objs as go
            from plotly.subplots import make_subplots
            # self.figures=[]

            # Symmetrized plots
            figure1 = make_subplots(rows=1, cols=1, subplot_titles=(['TMR']))
            figure2 = make_subplots(rows=1, cols=1, subplot_titles=(['MR']))
            figure3 = make_subplots(rows=1, cols=1, subplot_titles=(['Hall']))
            for data in self.symmetrized_data:
                color = int(
                    255.0
                    / len(self.symmetrized_data)
                    * self.symmetrized_data.index(data)
                )
                if data.rho_xx_up is not None and data.rho_xx_down is not None:
                    resistivity_tmr_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate((data.rho_xx_up, np.flip(data.rho_xx_down))),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure1.add_trace(resistivity_tmr_up, row=1, col=1)
                if data.mr_up is not None and data.mr_down is not None:
                    resistivity_mr_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate((data.mr_up, np.flip(data.mr_down))),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure2.add_trace(resistivity_mr_up, row=1, col=1)
                if data.rho_xy_up is not None:
                    resistivity_hall_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate((data.rho_xy_up, np.flip(data.rho_xy_down))),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure3.add_trace(resistivity_hall_up, row=1, col=1)
            figure1.update_layout(height=400, width=716, showlegend=True)
            figure2.update_layout(height=400, width=716, showlegend=True)
            figure3.update_layout(height=400, width=716, showlegend=True)
            self.figures.append(
                PlotlyFigure(label='TMR', figure=figure1.to_plotly_json())
            )
            self.figures.append(
                PlotlyFigure(label='MR', figure=figure2.to_plotly_json())
            )
            self.figures.append(
                PlotlyFigure(label='Hall', figure=figure3.to_plotly_json())
            )

            # Analyzed plots
            figure1 = make_subplots(rows=1, cols=1, subplot_titles=(['OHR']))
            figure2 = make_subplots(rows=1, cols=1, subplot_titles=(['AHR']))
            figure3 = make_subplots(rows=1, cols=1, subplot_titles=(['AHC']))
            for data in self.analyzed_data:
                color = int(
                    255.0 / len(self.analyzed_data) * self.analyzed_data.index(data)
                )
                if data.rho_ohe_up is not None and data.rho_ohe_down is not None:
                    resistivity_tmr_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate((data.rho_ohe_up, np.flip(data.rho_ohe_down))),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure1.add_trace(resistivity_tmr_up, row=1, col=1)
                if data.rho_ahe_up is not None and data.rho_ahe_down is not None:
                    resistivity_mr_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate((data.rho_ahe_up, np.flip(data.rho_ahe_down))),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure2.add_trace(resistivity_mr_up, row=1, col=1)
                if data.sigma_ahe_up is not None and data.sigma_ahe_down is not None:
                    resistivity_hall_up = go.Scatter(
                        x=np.concatenate((data.field, np.flip(data.field))),
                        y=np.concatenate(
                            (data.sigma_ahe_up, np.flip(data.sigma_ahe_down))
                        ),
                        name=data.title.split('at')[1].strip('.'),
                        marker_color=f'rgb({color},0,255)',
                        showlegend=True,
                    )
                    figure3.add_trace(resistivity_hall_up, row=1, col=1)
            figure1.update_layout(height=400, width=716, showlegend=True)
            figure2.update_layout(height=400, width=716, showlegend=True)
            figure3.update_layout(height=400, width=716, showlegend=True)
            self.figures.append(
                PlotlyFigure(label='OHR', figure=figure1.to_plotly_json())
            )
            self.figures.append(
                PlotlyFigure(label='AHR', figure=figure2.to_plotly_json())
            )
            self.figures.append(
                PlotlyFigure(label='AHC', figure=figure3.to_plotly_json())
            )


class CPFSPPMSACTMeasurement(CPFSPPMSMeasurement, PlotSection, EntryData):
    temperature_tolerance = Quantity(
        type=float,
        unit='kelvin',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='kelvin'
        ),
    )

    field_tolerance = Quantity(
        type=float,
        unit='gauss',
        a_eln=ELNAnnotation(component='NumberEditQuantity', defaultDisplayUnit='gauss'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        super().normalize(archive, logger)

        ### Start of the PPMSMeasurement normalizer

        if archive.data.sequence_file:
            logger.info('Parsing PPMS sequence file.')
            with archive.m_context.raw_file(self.sequence_file, 'r') as file:
                sequence = file.readlines()
                self.steps = find_ppms_steps_from_sequence(sequence)

        # Now create the according plots
        import plotly.express as px
        from plotly.subplots import make_subplots

        for data in self.data:
            if data.measurement_type == 'field':
                resistivity_ch1 = px.scatter(
                    x=data.magnetic_field, y=data.channels[0].resistivity
                )
                resistivity_ch2 = px.scatter(
                    x=data.magnetic_field, y=data.channels[1].resistivity
                )
            if data.measurement_type == 'temperature':
                resistivity_ch1 = px.scatter(
                    x=data.temperature, y=data.channels[0].resistivity
                )
                resistivity_ch2 = px.scatter(
                    x=data.temperature, y=data.channels[1].resistivity
                )
            figure1 = make_subplots(rows=2, cols=1, shared_xaxes=True)
            figure1.add_trace(resistivity_ch1.data[0], row=1, col=1)
            figure1.add_trace(resistivity_ch2.data[0], row=2, col=1)
            figure1.update_layout(height=400, width=716, title_text=data.name)
            self.figures.append(
                PlotlyFigure(label=data.name, figure=figure1.to_plotly_json())
            )


class CPFSPPMSMPMSMeasurement(CPFSPPMSMeasurement, PlotSection, EntryData):
    temperature_tolerance = Quantity(
        type=float,
        unit='kelvin',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='kelvin'
        ),
    )

    field_tolerance = Quantity(
        type=float,
        unit='gauss',
        a_eln=ELNAnnotation(component='NumberEditQuantity', defaultDisplayUnit='gauss'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        super().normalize(archive, logger)

        ### Start of the PPMSMeasurement normalizer

        if archive.data.sequence_file:
            logger.info('Parsing PPMS sequence file.')
            with archive.m_context.raw_file(self.sequence_file, 'r') as file:
                sequence = file.readlines()
                self.steps = find_ppms_steps_from_sequence(sequence)

        # Now create the according plots
        import plotly.express as px
        from plotly.subplots import make_subplots

        for data in self.data:
            if data.measurement_type == 'field':
                magnetization = px.scatter(x=data.magnetic_field, y=data.moment)
            if data.measurement_type == 'temperature':
                magnetization = px.scatter(x=data.temperature, y=data.moment)
            figure1 = make_subplots(rows=1, cols=1, shared_xaxes=True)
            figure1.add_trace(magnetization.data[0], row=1, col=1)
            figure1.update_layout(height=400, width=716, title_text=data.name)
            self.figures.append(
                PlotlyFigure(label=data.name, figure=figure1.to_plotly_json())
            )


class CPFSPPMSACMSMeasurement(CPFSPPMSMeasurement, PlotSection, EntryData):
    temperature_tolerance = Quantity(
        type=float,
        unit='kelvin',
        a_eln=ELNAnnotation(
            component='NumberEditQuantity', defaultDisplayUnit='kelvin'
        ),
    )

    field_tolerance = Quantity(
        type=float,
        unit='gauss',
        a_eln=ELNAnnotation(component='NumberEditQuantity', defaultDisplayUnit='gauss'),
    )

    def normalize(self, archive, logger: BoundLogger) -> None:  # noqa: PLR0912, PLR0915
        super().normalize(archive, logger)

        ### Start of the PPMSMeasurement normalizer

        if archive.data.sequence_file:
            logger.info('Parsing PPMS sequence file.')
            with archive.m_context.raw_file(self.sequence_file, 'r') as file:
                sequence = file.readlines()
                self.steps = find_ppms_steps_from_sequence(sequence)

        # Now create the according plots
        import plotly.express as px
        from plotly.subplots import make_subplots

        for data in self.data:
            if data.measurement_type == 'field':
                magnetization = px.scatter(x=data.magnetic_field, y=data.moment)
            if data.measurement_type == 'temperature':
                magnetization = px.scatter(x=data.temperature, y=data.moment)
            figure1 = make_subplots(rows=1, cols=1, shared_xaxes=True)
            figure1.add_trace(magnetization.data[0], row=1, col=1)
            figure1.update_layout(height=400, width=716, title_text=data.name)
            self.figures.append(
                PlotlyFigure(label=data.name, figure=figure1.to_plotly_json())
            )
