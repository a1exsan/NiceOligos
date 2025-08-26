from nicegui import app, ui
from raw_mat_widget import rawMatWidget
from raw_mat_widget import infoPanel
from raw_mat_widget import event

class rawmaterial_panel_page_model():
    def __init__(self, ip_addr, port):
        self.ip_addr, self.port = ip_addr, port
        self.pincode = ''

        self.set_widget_list()
        self.init_frontend()
        self.init_data()

    def set_widget_list(self):
        self.wi_list = {}
        self.wi_list['Solvents and reagents:'] = []
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000001', 'Ацетонитрил'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000007', 'Бензол'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000004', 'Пиридин'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000010', 'н-метилимидазол'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000018', 'проп. ангидрид'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000009', 'ДХУ'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000008', 'ДХУ'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000006', 'Ацетон'])
        self.wi_list['Solvents and reagents:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000011', 'Йод'])
        self.wi_list['Phosphoramidites:'] = []
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000120', 'dA амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000121', 'dC амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000122', 'dG амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000123', 'dT амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000129', 'FAM амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000133', 'Alkine амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000130', 'Biotin амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000124', 'lna A амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000126', 'lna C амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000125', 'lna G амидит'])
        self.wi_list['Phosphoramidites:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000127', 'lna T амидит'])
        self.wi_list['Azides:'] = []
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000167', 'R6G азид'])
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000168', 'SIMA азид'])
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000170', 'Cy5 азид'])
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000171', 'Cy5.5 азид'])
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000172', 'VIC азид'])

    def widgets_stack(self):
        wi_w = 200
        self.wi_list['widgets'] = {}
        with ui.row().classes("w-full"):
            self.stack_group_solvents = ui.label('Solvents and reagents:').style('font-size: 30px;')
            self.write_off_group_solvents = ui.button('Write-off', color='red', on_click=self.on_write_off_solvents)
            self.write_in_group_solvents = ui.button('Write-in', color='green', on_click=self.on_write_in_solvents)

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            for data in self.wi_list['Solvents and reagents:']:
                self.wi_list['widgets'][data[0]] = rawMatWidget(self.ip_addr, self.port, data[0], self.pincode, lbl=data[1])

        with ui.row().classes("w-full"):
            self.stack_group_amidits = ui.label('Phosphoramidites:').style('font-size: 30px;')
            self.write_off_group_amidits = ui.button('Write-off', color='red', on_click=self.on_write_off_amidits)
            self.write_in_group_amidits = ui.button('Write-in', color='green', on_click=self.on_write_in_amidits)

        wi_w = 200

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            for data in self.wi_list['Phosphoramidites:']:
                self.wi_list['widgets'][data[0]] = rawMatWidget(self.ip_addr, self.port, data[0], self.pincode, lbl=data[1])

        with ui.row().classes("w-full"):
            self.stack_group_azids = ui.label('Azides:').style('font-size: 30px;')
            self.write_off_group_azids = ui.button('Write-off', color='red', on_click=self.on_write_off_azids)
            self.write_in_group_azids = ui.button('Write-in', color='green', on_click=self.on_write_in_azids)

        wi_w = 200

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            for data in self.wi_list['Azides:']:
                self.wi_list['widgets'][data[0]] = rawMatWidget(self.ip_addr, self.port, data[0], self.pincode, lbl=data[1])

    def write_off_group(self, data, write_off=True):
        out = {}
        index = 1
        for row in data:
            out[index] = {1: row[1], 2: row[0]}
            index += 1
        self.info_panel.rawMat.selected_list = out.copy()
        if write_off:
            self.info_panel.on_write_off_stock(event())
        else:
            self.info_panel.on_write_in_stock(event())

    def on_write_off_solvents(self):
        self.write_off_group(self.wi_list['Solvents and reagents:'])

    def on_write_off_amidits(self):
        self.write_off_group(self.wi_list['Phosphoramidites:'])

    def on_write_off_azids(self):
        self.write_off_group(self.wi_list['Azides:'])

    def on_write_in_solvents(self):
        self.write_off_group(self.wi_list['Solvents and reagents:'], write_off=False)

    def on_write_in_amidits(self):
        self.write_off_group(self.wi_list['Phosphoramidites:'], write_off=False)

    def on_write_in_azids(self):
        self.write_off_group(self.wi_list['Azides:'], write_off=False)

    def init_frontend(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        if self.pincode != '':
            with (((ui.grid(columns=2).classes("w-full").style(f"grid-template-columns: {1550}px {1200}px")))):
                with ui.column():
                    self.widgets_stack()
                with ui.column():
                    self.lbl_info_panel = ui.label('Info panel:').style('font-size: 30px;')
                    self.info_panel = infoPanel(self.ip_addr, self.port, self.pincode, self.lbl_info_panel)
                    with ui.row():
                        ui.label('Search by name:').style('font-size: 30px;')
                        self.search_text = ui.input(label='', placeholder='Search name',
                                                    on_change=self.info_panel.rawMat.search_by_name).style('font-size: 30px;').classes('w-[700px]')
                        #self.on_search = ui.button('Search', color="#00a100").classes('w-[200px]')
                    #self.search_text.on_value_change = self.info_panel.search_by_name

        for data in self.wi_list['Solvents and reagents:']:
            self.wi_list['widgets'][data[0]].on_mouse_click = self.info_panel.on_show_info_menu_widget_click
        for data in self.wi_list['Phosphoramidites:']:
            self.wi_list['widgets'][data[0]].on_mouse_click = self.info_panel.on_show_info_menu_widget_click
        for data in self.wi_list['Azides:']:
            self.wi_list['widgets'][data[0]].on_mouse_click = self.info_panel.on_show_info_menu_widget_click

    def init_data(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')






