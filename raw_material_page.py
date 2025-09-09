from nicegui import app, ui
from raw_mat_widget import rawMatWidget
from raw_mat_widget import infoPanel
from raw_mat_widget import event
from raw_mat_widget import show_stock_operations
from OligoMap_utils import api_db_interface
import requests
import json
from datetime import datetime

class rawmaterial_panel_page_model(api_db_interface):
    def __init__(self, db_IP, db_port):
        super().__init__(db_IP, db_port)
        self.pincode = ''
        self.init_data()

        self.widget_db = 'gui_object_content_1.db'
        self.set_widget_list()
        self.init_frontend()


    def add_widgets_to_db(self, row):
        url = f"{self.api_db_url}/insert_data/{self.widget_db}/main_tab"
        r = requests.post(url,
                          json=json.dumps(
                              [
                                  row['object_id'], row['date'], row['date_format'], row['obj_json']
                              ]
                          )
                          , headers=self.headers())

    def update_widgets(self):
        row = {}
        row['#'] = self.wis_meta['#']
        row['object_id'] = self.wis_meta['object_id']
        row['date'] = self.wis_meta['date']
        row['date_format'] = self.wis_meta['date_format']
        row['obj_json'] = json.dumps(self.wi_list)

        url = f"{self.api_db_url}/update_data/{self.widget_db}/total_tab/{row['#']}"
        ret = requests.put(url, json=json.dumps(
                {
                    'name_list': ['object_id', 'date', 'date_format', 'obj_json'],
                    'value_list': [
                                    row['object_id'], row['date'], row['date_format'], row['obj_json']
                    ]
                }
            ), headers=self.headers())

    def load_widgets_content(self, object_id):
        url = f"{self.api_db_url}/get_keys_data/{self.widget_db}/main_tab/object_id/{object_id}"
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            return ret.json()
        else:
            return []

    def set_widget_list(self):
        self.wi_list = {}
        self.wis_meta = {}
        wis = self.load_widgets_content('raw_material_widgets')
        if len(wis) > 0:
            self.wis_meta['#'] = wis[0][0]
            self.wis_meta['object_id'] = wis[0][1]
            self.wis_meta['date'] = wis[0][2]
            self.wis_meta['date_format'] = wis[0][3]
            self.wi_list = json.loads(wis[0][4])
        else:
            self.init_widgets()

        #data = {}
        #data['object_id'] = 'raw_material_widgets'
        #data['date'] = datetime.now().date().strftime('%d.%m.%Y')
        #data['date_format'] = '%d.%m.%Y'
        #data['obj_json'] = json.dumps(self.wi_list)
        #print(data)
        #self.add_widgets_to_db(data)


    def init_widgets(self):
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
        self.wi_list['Azides:'].append(
            ['0ce92bd9f0d44e501389eff9f5a803d6ddd8a650fffc58043f633f9183cf91dd', 'TAMRA азид'])
        self.wi_list['Azides:'].append(['b2db68279078604bc7783320700cd2dccb3c16c72798a3e72ac5e51568958cab', 'Cy7 азид'])
        self.wi_list['Azides:'].append(['INIT_BASE_CODE_OLIGO_LAB_0000166', 'ROX азид'])


    def widgets_stack(self):
        self.wi_key_list = list(self.wi_list.keys())
        self.wi_list['widgets'] = {}
        for key in self.wi_key_list:
            with ui.row().classes("w-full"):
                ui.label(key).style('font-size: 30px;')
                ui.button(f'Write-off_{key}', color='red', on_click=self.on_write_off_group)
                ui.button(f'Write-in_{key}', color='green', on_click=self.on_write_in_group)
                ui.button(f'Delete_{key}', color='orange', on_click=self.on_delete_wi_group)

            wi_w = 200

            with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
                for data in self.wi_list[key]:
                    self.wi_list['widgets'][data[0]] = rawMatWidget(self.db_IP, self.db_port, data[0], self.pincode, lbl=data[1])


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

    def delete_wi_group(self, group, data):
        out = {}
        index = 1
        for row in data:
            out[index] = {1: row[1], 2: row[0]}
            index += 1
        self.info_panel.rawMat.selected_list = out.copy()
        self.info_panel.on_delete_wi_in_stock(
            {
                'group': group,
                'group_data': data,
                'wi_list': self.wi_list,
                'meta': self.wis_meta
            }
        )


    def on_write_off_group(self, e):
        s = e.sender.props['label']
        wis_group = s[s.find('_')+1:]
        if wis_group in list(self.wi_list.keys()):
            self.write_off_group(self.wi_list[wis_group])

    def on_write_in_group(self, e):
        s = e.sender.props['label']
        wis_group = s[s.find('_') + 1:]
        if wis_group in list(self.wi_list.keys()):
            self.write_off_group(self.wi_list[wis_group], write_off=False)

    def on_delete_wi_group(self, e):
        s = e.sender.props['label']
        wis_group = s[s.find('_') + 1:]
        if wis_group in list(self.wi_list.keys()):
            self.delete_wi_group(wis_group, self.wi_list[wis_group])

    def init_frontend(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        if self.pincode != '':
            with (((ui.grid(columns=2).classes("w-full").style(f"grid-template-columns: {1550}px {1200}px")))):
                with ui.column():
                    self.widgets_stack()
                with ui.column():
                    self.lbl_info_panel = ui.label('Info panel:').style('font-size: 30px;')
                    self.info_panel = infoPanel(self.db_IP, self.db_port, self.pincode, self.lbl_info_panel)
                    with ui.row():
                        ui.label('Search by name:').style('font-size: 30px;')
                        self.search_text = ui.input(label='', placeholder='Search name',
                                                    on_change=self.info_panel.rawMat.search_by_name).style('font-size: 30px;').classes('w-[700px]')
                    ui.button(f'show operations', color='green', on_click=self.on_show_write_of_in_operations)
                        #self.on_search = ui.button('Search', color="#00a100").classes('w-[200px]')
                    #self.search_text.on_value_change = self.info_panel.search_by_name

            for key in self.wi_key_list:
                for data in self.wi_list[key]:
                    self.wi_list['widgets'][data[0]].on_mouse_click = self.info_panel.on_show_info_menu_widget_click

    def init_data(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')

    def on_show_write_of_in_operations(self):
        stock_dialog = show_stock_operations(self.db_IP, self.db_port)
        stock_dialog.dialog.open()






