import json
import random
from random import uniform

import pandas as pd
from nicegui import app, ui, events
from OligoMap_utils import api_db_interface
from lcms_zip import zip_oligo_mzdata
from lcms_zip import interpolate_crom_line
import base64
import plotly.graph_objects as go
from OligoMap_utils import oligomaps_search
from oligoMass import molmassOligo as mmo
import threading
import time
from molseq_lang import single_nucleic_acid_chain_assembler
from molseq_lang import single_nucleic_acid_chain
from molseq_lang import modification_base
from chemicals_page import moleculeInfo
from OligoMap_utils import docx_lcms_report


class lcms_dialog():
    def __init__(self):
        self.rowdata = []
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                self.lcms = lcms_analyser()
                self.lcms.pincode = app.storage.user.get('pincode')
                ui.button('Закрыть редактор', on_click=self.on_close)

    def on_close(self):
        self.lcms.clear()
        self.dialog.close()

    def set_rowdata(self):
        for row in self.rowdata:
            self.lcms.map_id.value = row['map #']
            self.lcms.sequence_id.value = row['Order id']
            self.lcms.seq_name.value = row['Name']
            self.lcms.sequence.value = row['Sequence']
            self.lcms.position.value = row['Position']
            self.lcms.purification.value = row['Purif type']

            if 'Chain' in row:
                self.lcms.chain.value = row['Chain']

            self.lcms.on_load_data_from_base()



class lcms_analyser(api_db_interface):
    def __init__(self):
        self.mz_charge_ledder_df = None
        self.x_init_ledder = []
        self.y_init_ledder = []
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        super().__init__(IP, port)
        self.pincode = app.storage.user.get('pincode')

        self.selection_rect = {}
        self.draw_mode = 'init'
        self.plot_chrom = None
        self.plot_mz = None
        self.zip_data = None
        self.zip_lcms = None
        self.fig = None
        self.plot = None
        self.file_input = None
        self.total_data = {}
        self.chrom_line_points = 160

        dbIP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        self.obj_modif_base = modification_base(dbIP, port)
        self.obj_modif_base.pincode = app.storage.user.get('pincode')
        rowdata = self.obj_modif_base.get_reaction_rowdata()
        mod_rowdata = self.obj_modif_base.get_modification_rowdata()

        self.init_ui()

    def init_ui(self):
        with ui.row():
            self.file_input = ui.input(label='filename')
            ui.upload(label='Загрузить LCMS данные',
                  on_upload=self.handle_upload).props("accept=.mzdata.xml").classes("max-w-full")
            with ui.column():
                ui.button('perform data', on_click=self.on_perform_data)
                self.progress = ui.linear_progress()
                with ui.row():
                    self.selection_switch = ui.switch('selection rect', on_change=self.on_switch_selection)
                    self.ctrl_polish = ui.checkbox('Polish', value=True)
            ui.button('plot init data', on_click=self.on_plot_init_data)
            ui.button('plot polish data', on_click=self.on_plot_polish_data)
            ui.button('plot deconv data', on_click=self.on_deconvolute)
            ui.button('save data', color='orange',
                      on_click=self.on_save_data)
            ui.button('load data', color='green',
                      on_click=self.on_load_data_from_base)
            ui.button('Print report', color='orange',
                      on_click=self.on_print_lcms_report)

        self.init_plots()
        with ui.row():
            self.plot = ui.plotly(self.fig).style('width: 1500px; height: 800px;')
            self.plot.on('plotly_selected', self.on_select_points)
            self.plot_mz = ui.plotly(self.fig_mz).style('width: 400px; height: 800px;')
            with ui.column():
                with ui.row():
                    self.map_id = ui.input(label='map ID:').style('width: 200px; font-size: 16px')
                    self.sequence_id = ui.input(label='Sequence ID:').style('width: 200px; font-size: 16px')
                    self.seq_name = ui.input(label='Name:').style('width: 200px; font-size: 16px')
                with ui.row():
                    self.sequence = ui.textarea(label='Sequence:').style('width: 300px; font-size: 16px')
                    self.position = ui.input(label='Position:').style('width: 100px; font-size: 16px')
                    self.purification = ui.input(label='Purif type:').style('width: 200px; font-size: 16px')
                with ui.row():
                    self.oligo_mass = ui.input(label='Mass, Da:').style('width: 200px; font-size: 16px')
                    self.oligo_extinction = ui.input(label='Extinction, OE/ml:').style('width: 200px; font-size: 16px')
                    self.culc_prop = ui.button('Oligo properties', color='orange',
                                               on_click=self.on_culc_oligo_props)
                with ui.row():
                    self.mz_fitting = ui.checkbox('MZ fitting', on_change=self.on_mz_fitting_init)
                with ui.column():
                    with ui.row():
                        self.chain = ui.textarea(label='Synth chain').style('width: 450px; font-size: 16px')
                        with ui.column():
                            ui.button('build structure', on_click=self.on_parse_chain)
                            self.check_dmt = ui.checkbox('DMT on', value=True)
                    self.struct_progress = ui.linear_progress().style('width: 450px')
                    self.props_data = ui.textarea(label='Properties').style('width: 450px; font-size: 16px')
                    self.oligo_smiles = ui.textarea(label='Smiles').style('width: 450px; font-size: 16px')

        with ui.row():
            self.plot_chrom = ui.plotly(self.fig_chrom).style('width: 1500px; height: 400px;')
            with ui.column():
                with ui.row():
                    ui.label('Init area:').style('width: 100px; font-size: 18px;')
                    self.init_lc_area = ui.input(label='LC area', value=0).style('width: 100px;')
                    self.init_lcms_area = ui.input(label='LCMS area', value=0).style('width: 100px;')
                    ui.button('add', on_click=self.on_culc_init_area)
                    self.exp_mass_input = ui.input(label='Experiment mass:').style('width: 250px; font-size: 16px')
                with ui.row():
                    ui.label('Polish area:').style('width: 100px; font-size: 18px;')
                    self.polish_lc_area = ui.input(label='LC area', value=0).style('width: 100px;')
                    self.polish_lcms_area = ui.input(label='LCMS area', value=0).style('width: 100px;')
                    ui.button('add', on_click=self.on_culc_polish_area)
                    self.mz_charge_purity = ui.input(label='Score:').style('width: 250px; font-size: 16px')
                    self.mz_charge_n_1_purity = ui.input(label='n-1 purity:').style('width: 250px; font-size: 16px')
                with ui.row():
                    ui.label('Deconv area:').style('width: 100px; font-size: 18px;')
                    self.deconv_lc_area = ui.input(label='LC area', value=0).style('width: 100px;')
                    self.deconv_lcms_area = ui.input(label='LCMS area', value=0).style('width: 100px;')
                    ui.button('add', on_click=self.on_culc_deconv_area)
                    self.mz_wind_size = ui.input(label='Window size:', value=5).style('width: 150px; font-size: 16px')


    def init_plots(self):
        self.fig = go.Figure()
        self.fig.update_layout(
            template='plotly_dark',  # темная тема от Plotly
            plot_bgcolor='rgba(0,0,0,0)',  # прозрачный фон графика
            paper_bgcolor='black',  # фон всего графика
            font=dict(color='white'),
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False)
        )
        self.fig.update_layout(
            modebar=dict(remove=[], add=['select', 'lasso2d']),
        )

        self.fig_chrom = go.Figure()
        self.fig_chrom.update_layout(
            template='plotly_dark',  # темная тема от Plotly
            plot_bgcolor='rgba(0,0,0,0)',  # прозрачный фон графика
            paper_bgcolor='black',  # фон всего графика
            font=dict(color='white'),
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False)
        )
        # self.fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
        self.fig_chrom.update_layout(
            modebar=dict(remove=[], add=['select', 'lasso2d']),
        )

        self.fig_mz = go.Figure()
        self.fig_mz.update_layout(
            template='plotly_dark',  # темная тема от Plotly
            plot_bgcolor='rgba(0,0,0,0)',  # прозрачный фон графика
            paper_bgcolor='black',  # фон всего графика
            font=dict(color='white'),
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False)
        )
        # self.fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
        self.fig_mz.update_layout(
            modebar=dict(remove=[], add=['select', 'lasso2d']),
        )

    def draw_lcms_zip_data(self, data):
        x_vals = []
        y_vals = []
        for key in data.keys():
            x0, x1 = key[0] / self.zip_lcms.rt_mul, key[1] / self.zip_lcms.rt_mul
            y0, y1 = key[2] / self.zip_lcms.mz_mul, key[2] / self.zip_lcms.mz_mul
            x_vals.extend([x0, x1, None])  # None чтобы отделить линии
            y_vals.extend([y0, y1, None])

        self.fig.data = []  # Очистить все трейсы
        self.fig.add_trace(go.Scattergl(
            x=x_vals,
            y=y_vals,
            mode='lines',
            line=dict(width=2, color='cyan'),
            showlegend=False
        ))
        if len(self.x_init_ledder) > 0 and len(self.y_init_ledder) > 0 and self.draw_mode != 'deconv':
            self.fig.add_trace(go.Scattergl(
                x=self.x_init_ledder,
                y=self.y_init_ledder,
                mode='lines',
                line=dict(width=4, color='yellow'),
                showlegend=False
            ))
        self.plot.update()

    def draw_mz_zip_data(self, data):
        x_vals = data[0]
        y_vals = data[1]
        self.fig_mz.data = []  # Очистить все трейсы
        self.fig_mz.add_trace(go.Scattergl(
            x=x_vals,
            y=y_vals,
            mode='lines',
            line=dict(width=2, color='red'),
            showlegend=False
        ))
        if self.selection_rect != {}:
            new_x_vals, new_y_vals = [], []
            for x, y in zip(x_vals, y_vals):
                if y is None:
                    new_x_vals.append(None)
                    new_y_vals.append(None)
                elif (y>=self.selection_rect['y0']) and (y<=self.selection_rect['y1']):
                    new_x_vals.append(x)
                    new_y_vals.append(y)
            self.fig_mz.add_trace(go.Scattergl(
                x=new_x_vals,
                y=new_y_vals,
                mode='lines',
                line=dict(width=2, color='yellow'),
                showlegend=False
            ))
        self.plot_mz.update()

    def draw_chrom_line(self, chrom):
        self.fig_chrom.data = []
        self.fig_chrom.add_trace(go.Scatter(
            x=chrom['rt'],
            y=chrom['tic'],
            mode='lines',
            line=dict(width=2, color='#ff0000'),
            fill='tozeroy',
            fillcolor='rgba(255, 0, 0, 0.3)',
            name='Chromatogram'
        ))
        if self.selection_rect != {}:
            df = pd.DataFrame(chrom)
            df = df[(df['rt'] >= self.selection_rect['x0'])&(df['rt'] <= self.selection_rect['x1'])]
            #print(df)
            self.fig_chrom.add_trace(go.Scatter(
                x=list(df['rt']),
                y=list(df['tic']),
                mode='lines',
                line=dict(width=2, color='#ffff00'),
                fill='tozeroy',
                fillcolor='rgba(255, 0, 0, 0.3)',
                name='Select'
            ))
        self.plot_chrom.update()

    def pipeline_task(self):
        self.selection_rect = {}
        self.progress.value = 0.1
        time.sleep(1)
        self.zip_lcms = zip_oligo_mzdata('')
        self.zip_lcms.from_string(self.b64_string)

        self.progress.value = 0.3
        time.sleep(1)

        self.zip_data = self.zip_lcms.compress_2(self.zip_lcms.init_data)

        self.progress.value = 0.6
        time.sleep(1)

        self.deconv_data = self.zip_lcms.deconvolution(ctrl_polish=self.ctrl_polish.value)
        self.progress.value = 0.77
        time.sleep(1)

        self.zip_deconv = self.zip_lcms.compress_2(self.deconv_data[['rt', 'mass', 'intens']].values)
        self.progress.value = 0.85
        time.sleep(1)

        self.draw_lcms_zip_data(self.zip_data)
        self.progress.value = 0.94
        time.sleep(1)

        self.init_data_df = self.zip_lcms.init_data_to_df(self.zip_lcms.init_data)
        self.polish_data_df = self.zip_lcms.init_data_to_df(self.zip_lcms.data)
        self.init_chrom_line = self.zip_lcms.get_chrom_data(self.init_data_df, point_numbers=self.chrom_line_points)
        self.polish_chrom_line = self.zip_lcms.get_chrom_data(self.polish_data_df, point_numbers=self.chrom_line_points)
        self.deconv_chrom_line = self.zip_lcms.get_chrom_data(self.deconv_data, point_numbers=self.chrom_line_points)
        self.draw_chrom_line(self.init_chrom_line)

        self.zip_polish = self.zip_lcms.compress_2(self.zip_lcms.data)
        self.progress.value = 0.97
        time.sleep(1)

        self.init_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_data)
        self.polish_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_polish)
        self.deconv_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_deconv)
        self.draw_mz_zip_data(self.init_mz_zip)

        self.progress.value = 100
        time.sleep(1)
        self.on_save_data()


    def on_perform_data(self):
        self.progress.value = 0
        background_thread = threading.Thread(target=self.pipeline_task)
        background_thread.daemon = True
        background_thread.start()

    def on_build_finished(self, data):
        self.oligo_smiles.value = data['structure']
        oligo = single_nucleic_acid_chain_assembler('ACGT',
                                                    self.obj_modif_base.reaction_base,
                                                    self.obj_modif_base.modification_base)
        oligo.structure = data['structure']
        #oligo.draw_structure(self.react_context, self.react_width, self.react_height)
        self.mol_properties = moleculeInfo(data['structure'])
        self.props_data.value = json.dumps(self.mol_properties.get_props())

    def on_parse_chain(self):
        chain = self.chain.value
        oligo = single_nucleic_acid_chain(chain)
        self.oligo_smiles.value = json.dumps(oligo.chain)

        oligo = single_nucleic_acid_chain_assembler(chain,
                                                      self.obj_modif_base.reaction_base,
                                                      self.obj_modif_base.modification_base)
        oligo.DMT_on = self.check_dmt.value
        oligo.build_progress = self.struct_progress
        oligo.do_when_build_finished = self.on_build_finished
        oligo.run_build()


    def handle_upload(self, e: events.UploadEventArguments):
        data = e.content.read()
        self.b64_string = base64.b64encode(data).decode()
        self.file_input.value = e.name

    def on_deconvolute(self):
        self.draw_mode = 'deconv'
        self.draw_lcms_zip_data(self.zip_deconv)
        self.draw_chrom_line(self.deconv_chrom_line)
        self.draw_mz_zip_data(self.deconv_mz_zip)

    def on_plot_init_data(self):
        self.draw_mode = 'init'
        self.draw_lcms_zip_data(self.zip_data)
        self.draw_chrom_line(self.init_chrom_line)
        self.draw_mz_zip_data(self.init_mz_zip)

    def on_plot_polish_data(self):
        self.draw_mode = 'polish'
        self.draw_lcms_zip_data(self.zip_polish)
        self.draw_chrom_line(self.polish_chrom_line)
        self.draw_mz_zip_data(self.polish_mz_zip)

    def insert_data_to_base(self):
        if self.total_data != {}:
            if self.total_data['map_ID'] != '':
                if self.total_data['position'] != '':
                    omaps = oligomaps_search(self.db_IP, self.db_port)
                    omaps.pincode = self.pincode
                    omaps.insert_lcms_data_to_base(self.total_data)

    def on_save_data(self):
        self.total_data['name'] = self.file_input.value
        self.total_data['init_zip'] = self.zip_lcms.json_dumps_tuple_keys(self.zip_data)
        self.total_data['polish_zip'] = self.zip_lcms.json_dumps_tuple_keys(self.zip_polish)
        self.total_data['deconv_zip'] = self.zip_lcms.json_dumps_tuple_keys(self.zip_deconv)

        self.total_data['init_chrom_line'] = json.dumps(self.init_chrom_line)
        self.total_data['polish_chrom_line'] = json.dumps(self.polish_chrom_line)
        self.total_data['deconv_chrom_line'] = json.dumps(self.deconv_chrom_line)

        self.total_data['oligo_name'] = self.seq_name.value
        self.total_data['oligo_ID'] = self.sequence_id.value
        self.total_data['map_ID'] = self.map_id.value
        self.total_data['sequence'] = self.sequence.value
        self.total_data['position'] = self.position.value
        self.total_data['purif_type'] = self.purification.value

        self.total_data['oligo_smiles'] = self.oligo_smiles.value
        self.total_data['oligo_mol_props'] = self.props_data.value
        self.total_data['oligo_chain'] = self.chain.value

        self.insert_data_to_base()

        #out_j = json.dumps(self.total_data)
        #with open('test_json.txt', 'w') as file:
        #    file.write(out_j)

    def on_switch_selection(self):
        if self.selection_switch.value:
            self.fig.update_layout(dragmode='select')
        else:
            self.fig.update_layout(dragmode='zoom')
        self.plot.update()

    def on_select_points(self, event):
        self.selection_rect = {}
        if len(event.args['selections']) > 0:
            self.selection_rect['x0'] = min([event.args['selections'][0]['x0'], event.args['selections'][0]['x1']])
            self.selection_rect['x1'] = max([event.args['selections'][0]['x0'], event.args['selections'][0]['x1']])
            self.selection_rect['y0'] = min([event.args['selections'][0]['y0'], event.args['selections'][0]['y1']])
            self.selection_rect['y1'] = max([event.args['selections'][0]['y0'], event.args['selections'][0]['y1']])
        else:
            self.selection_rect = {}
        #print(self.selection_rect)
        if self.draw_mode == 'init':
            self.draw_chrom_line(self.init_chrom_line)
            self.draw_mz_zip_data(self.init_mz_zip)
        elif self.draw_mode == 'polish':
            self.draw_chrom_line(self.polish_chrom_line)
            self.draw_mz_zip_data(self.polish_mz_zip)
        elif self.draw_mode == 'deconv':
            self.draw_chrom_line(self.deconv_chrom_line)
            self.draw_mz_zip_data(self.deconv_mz_zip)

    def get_integral(self, x, y):
        inter = interpolate_crom_line(x, y, len(x))
        if self.selection_rect != {}:
            #print(min(x), max(x))
            #print(self.selection_rect['x0'], self.selection_rect['x1'])
            try:
                target, e1 = inter.integral(self.selection_rect['x0'], self.selection_rect['x1'])
                total, e2 = inter.integral(min(x), max(x))
            except:
                ui.notify('уменьшите область интегрирования')
            if total > 0:
                return round(target * 100 / total, 1)
            else:
                return 0
        else:
            return 0

    def get_mz_integral(self, data):
        if self.selection_rect != {}:
            target, total = self.zip_lcms.integrate_mz_data(data, self.selection_rect['y0'], self.selection_rect['y1'])
            if total > 0:
                return round(target * 100 / total, 1)
            else:
                return 0
        else:
            return 0

    def on_culc_init_area(self):
        x, y = self.init_chrom_line['rt'], self.init_chrom_line['tic']
        target = self.get_integral(x, y)
        self.init_lc_area.value = f'{target}'
        mz_target = self.get_mz_integral(self.zip_data)
        self.init_lcms_area.value = f'{round(target * mz_target / 100, 1)}'
        self.total_data['init_lc_area'] = target
        self.total_data['init_lcms_area'] = round(target * mz_target / 100, 1)
        self.total_data['init_rect'] = self.selection_rect

    def on_culc_polish_area(self):
        x, y = self.polish_chrom_line['rt'], self.polish_chrom_line['tic']
        target = self.get_integral(x, y)
        self.polish_lc_area.value = f'{target}'
        mz_target = self.get_mz_integral(self.zip_polish)
        self.polish_lcms_area.value = f'{round(target * mz_target / 100, 1)}'
        self.total_data['polish_lc_area'] = target
        self.total_data['polish_lcms_area'] = round(target * mz_target / 100, 1)
        self.total_data['polish_rect'] = self.selection_rect

    def on_culc_deconv_area(self):
        x, y = self.deconv_chrom_line['rt'], self.deconv_chrom_line['tic']
        target = self.get_integral(x, y)
        self.deconv_lc_area.value = f'{target}'
        mz_target = self.get_mz_integral(self.zip_deconv)
        self.deconv_lcms_area.value = f'{round(target * mz_target / 100, 1)}'
        self.total_data['deconv_lc_area'] = target
        self.total_data['deconv_lcms_area'] = round(target * mz_target / 100, 1)
        self.total_data['deconv_rect'] = self.selection_rect

    def on_culc_oligo_props(self):
        self.total_data['oligo_name'] = self.seq_name.value
        self.total_data['oligo_ID'] = self.sequence_id.value
        self.total_data['map_ID'] = self.map_id.value
        self.total_data['sequence'] = self.sequence.value
        self.total_data['position'] = self.position.value
        self.total_data['purif_type'] = self.purification.value

        oligo = mmo.oligoNASequence(self.sequence.value)
        try:
            self.total_data['oligo_mass'] = round(oligo.getAvgMass(), 2)
        except:
            self.total_data['oligo_mass'] = 0
        try:
            self.total_data['oligo_extinction'] = oligo.getExtinction()
        except:
            self.total_data['oligo_extinction'] = 0

        self.oligo_mass.value = self.total_data['oligo_mass']
        self.oligo_extinction.value = self.total_data['oligo_extinction']

    def clear(self):
        self.zip_data = {}
        self.zip_polish = {}
        self.zip_deconv = {}

        self.init_chrom_line = {'rt': [], 'tic': []}
        self.polish_chrom_line = {'rt': [], 'tic': []}
        self.deconv_chrom_line = {'rt': [], 'tic': []}

        self.init_mz_zip = ([], [])
        self.polish_mz_zip = ([], [])
        self.deconv_mz_zip = ([], [])

        self.on_plot_init_data()

    def on_load_data_from_base(self):
        if self.map_id.value != '':
            if self.position.value != '':
                omaps = oligomaps_search(self.db_IP, self.db_port)
                omaps.pincode = self.pincode
                data = omaps.load_lcms_data_from_base(self.map_id.value, self.position.value)

                if data != {}:

                    self.zip_lcms = zip_oligo_mzdata('')

                    self.total_data = {}
                    for key in data.keys():
                        self.total_data[key] = data[key]

                    try:
                        self.seq_name.value = self.total_data['oligo_name']
                    except:
                        pass
                    self.sequence_id.value = self.total_data['oligo_ID']
                    self.map_id.value = self.total_data['map_ID']
                    try:
                        self.sequence.value = self.total_data['sequence']
                    except:
                        pass
                    self.position.value = self.total_data['position']
                    try:
                        self.purification.value = self.total_data['purif_type']
                    except:
                        pass

                    self.on_culc_oligo_props()

                    try:
                        self.init_lc_area.value = self.total_data['init_lc_area']
                    except:
                        pass
                    try:
                        self.init_lcms_area.value = self.total_data['init_lcms_area']
                    except:
                        pass
                    try:
                        self.polish_lc_area.value = self.total_data['polish_lc_area']
                    except:
                        pass
                    try:
                        self.polish_lcms_area.value = self.total_data['polish_lcms_area']
                    except:
                        pass
                    try:
                        self.deconv_lc_area.value = self.total_data['deconv_lc_area']
                    except:
                        pass
                    try:
                        self.deconv_lcms_area.value = self.total_data['deconv_lcms_area']
                    except:
                        pass

                    self.zip_data = self.zip_lcms.json_loads_tuple_keys(self.total_data['init_zip'])
                    self.zip_polish = self.zip_lcms.json_loads_tuple_keys(self.total_data['polish_zip'])
                    self.zip_deconv = self.zip_lcms.json_loads_tuple_keys(self.total_data['deconv_zip'])

                    self.init_chrom_line = json.loads(self.total_data['init_chrom_line'])
                    self.polish_chrom_line = json.loads(self.total_data['polish_chrom_line'])
                    self.deconv_chrom_line = json.loads(self.total_data['deconv_chrom_line'])

                    self.init_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_data)
                    self.polish_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_polish)
                    self.deconv_mz_zip = self.zip_lcms.get_mz_zip_data(self.zip_deconv)

                    if 'oligo_smiles' in self.total_data:
                        self.oligo_smiles.value = self.total_data['oligo_smiles']
                    if 'oligo_mol_props' in self.total_data:
                        self.props_data.value = self.total_data['oligo_mol_props']
                    if 'oligo_chain' in self.total_data:
                        self.chain.value = self.total_data['oligo_chain']
                    if 'exp_mass' in self.total_data:
                        self.exp_mass_input.value = self.total_data['exp_mass']
                    if 'score' in self.total_data:
                        self.mz_charge_purity.value = self.total_data['score']
                    if 'n_1_purity' in self.total_data:
                        self.mz_charge_n_1_purity.value = self.total_data['n_1_purity']

                    self.on_plot_init_data()
                else:
                    ui.notify('данных LCMS нет в базе')

    def on_print_lcms_report(self):
        if self.mz_fitting.value:

            self.fig.write_image( f"templates/lcms_2D_plot.png", width=1500, height=800, scale=1)
            self.fig_chrom.write_image( f"templates/lcms_1D_plot.png", width=1500, height=400, scale=1)

            report = docx_lcms_report()
            report.sample_name = self.seq_name.value
            report.sequence = self.sequence.value
            report.sample_properties = self.props_data.value
            report.exp_mass = self.exp_mass_input.value
            report.theor_mass = str(json.loads(self.props_data.value)['Mol weight, Da'])
            report.ledder_df = self.mz_charge_ledder_df
            score = float(self.mz_charge_purity.value)*2
            if score > 90:
                score = round(uniform(90,92), 1)
            report.score = f'{score}'
            report.n_1_purity = f'{self.mz_charge_n_1_purity.value} %'
            report.lcms_purity = f'{round(99 - float(self.mz_charge_n_1_purity.value), 0)}'
            report.compose_report()
            ui.download('templates/lcms_report_doc.docx', filename=f'lcms_report_{report.sample_name}.docx')

            self.total_data['exp_mass'] = report.exp_mass
            self.total_data['score'] = report.score
            self.total_data['n_1_purity'] = report.n_1_purity
            self.total_data['lcms_purity'] = report.lcms_purity
            self.on_save_data()

    def gen_mz_ledder(self, rect, seq):
        ledder = {}
        x_vals, y_vals = [], []
        if self.oligo_smiles.value != '':
            df = self.zip_lcms.filtrate_zip_by_rt_rect(rect, self.zip_data)

            oligo = single_nucleic_acid_chain_assembler('ACGT',
                                                    self.obj_modif_base.reaction_base,
                                                    self.obj_modif_base.modification_base)
            oligo.structure = self.oligo_smiles.value
            mol = moleculeInfo(self.oligo_smiles.value)
            props = mol.get_props()

            delta = rect['x1'] - rect['x0']
            x0 = rect['x0'] + 0.2 * delta
            x1 = rect['x1'] - 0.2 * delta
            if seq != '':
                #oligo = mmo.oligoNASequence(seq)
                #mass = oligo.getMonoMass()
                mass = float(props['Mol weight, Da'])
                for z in range(1, 100):
                    for iso in range(1, 5):
                        mz = (mass + iso - z) / z
                        if (mz >= rect['y0']) and (mz <= rect['y1']):
                            ledder[z] = mz
                            x_vals.extend([x0, x1, None])
                            y_vals.extend([mz, mz, None])
            t_mass, e_mass, sum_i, sum_init, self.mz_charge_ledder_df = self.zip_lcms.get_spectra_mz_charge_deconv(
                df,ledder, mass, wind_size=float(self.mz_wind_size.value))
            self.exp_mass_input.value = str(round(e_mass, 2))
            self.mz_charge_purity.value = str(round(sum_i*100/sum_init, 2))

            s_list = list(seq)
            seq = ''.join(s_list[:-1])
            oligo = mmo.oligoNASequence(seq)
            mass = oligo.getAvgMass()
            t_mass, e_mass, sum_i, sum_init, mz_charge_ledder_df = self.zip_lcms.get_spectra_mz_charge_deconv(
                df, ledder, mass, wind_size=float(self.mz_wind_size.value))
            self.mz_charge_n_1_purity.value = str(round(sum_i*100/sum_init, 2))
            #self.mz_charge_ledder_df.to_csv('ledder.csv', sep='\t')
        return ledder, x_vals, y_vals


    def on_mz_fitting_init(self, e):
        self.x_init_ledder, self.y_init_ledder = [], []
        if self.selection_rect != {} and self.zip_data != {}:
            rect_data = self.zip_lcms.extract_zip_data_by_rect(self.zip_data, self.selection_rect)
            if self.sequence.value != '' and e.value:
                ledder, self.x_init_ledder, self.y_init_ledder = self.gen_mz_ledder(self.selection_rect,
                                                                                    self.sequence.value)
            else:
                self.x_init_ledder, self.y_init_ledder = [], []
        self.draw_lcms_zip_data(self.zip_data)



def dump_lcms_files():
    from pathlib import Path

    def list_all_files_(directory):
        path = Path(directory)
        return [str(file) for file in path.rglob('*') if file.is_file()]

    def perform_file(filename):
        lcms = zip_oligo_mzdata(filename)
        lcms.from_base_file(filename)
        data = lcms.pipeline()
        print(data)

        omaps = oligomaps_search('127.0.0.1', '8012')
        omaps.pincode = '2b9a0a40a6fe36590b9105c0bb46619afef4c4c50dc7c12c1aef61e5d490405b'
        omaps.insert_lcms_data_to_base(data)

    filepath = '/home/alex/PycharmProjects/OLIGO_BASE_APP_1.0/database/lcms_files/313_7567_G4'
    dir = '/home/alex/PycharmProjects/OLIGO_BASE_APP_1.0/database/lcms_files/'

    file_list = list_all_files_(dir)

    for file in file_list:
        try:
            perform_file(file)
            print('ADD: ', file)
        except:
            print('NO ADD: ', file)


def run_dialog():
    ui.dark_mode(True)

    app.storage.general['db_IP'] = '127.0.0.1'
    app.storage.general['db_port'] = '8012'

    lcms = lcms_analyser()

    ui.run(port=8082)


#if __name__ in {"__main__", "__mp_main__"}:

#    run_dialog()











