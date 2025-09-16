from nicegui import events, ui, app
from OligoMap_utils import api_db_interface, oligomaps_search
import requests
import pandas as pd
import json
import datetime
from chemicals_page import phys_chem_props_interface
import hashlib



class show_stock_operations(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        else:
            self.pincode = ''

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        tab_df = pd.DataFrame({
            '#': [],
            'Name': [],
            "Amount": [],
            'Unicode': [],
            "Date": [],
            "Time": [],
            'User': []
        })

        colDefs = [
            {"field": "#", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Name", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Amount", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Unicode", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Date", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Time", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "User", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
        ]

        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                self.ag_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': tab_df.to_dict('records'),
                        'rowSelection': 'multiple',
                        "pagination": True,
                        # "enableRangeSelection": True,
                    }
                    ,
                    theme='alpine-dark').style('width: 1200px; height: 800px')  # alpine  material  quartz  balham
                self.ag_grid.auto_size_columns = True

                with ui.row():
                    ui.button('Списания', on_click=self.on_write_off)
                    ui.button('Поступления', on_click=self.on_write_in)


    def set_ag_grid_tab(self, name):

        users = self.get_all_data_in_tab(f'users')
        ids = {}
        for user in users:
            ids[user[2]] = user[1]

        data = self.get_all_data_in_tab(f'{name}')
        df = pd.DataFrame(data)
        tab_df = pd.DataFrame({
            '#': df[0],
            'Name': df[1],
            "Amount": df[2],
            'Unicode': df[3],
            "Date": df[4],
            "Time": df[5],
            'User': df[6]
        })
        tab_df['User'] = [ids[i] for i in tab_df['User']]
        self.ag_grid.options['rowData'] = tab_df.to_dict('records')
        self.ag_grid.update()

    def on_write_off(self):
        self.set_ag_grid_tab('output_tab')

    def on_write_in(self):
        self.set_ag_grid_tab('input_tab')

    def get_all_data_in_tab(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/{tab_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()


class reagent_form_dialog():
    def __init__(self, data, new=False):
        self.data = data

        if new:
            self.data_obj = {}
            desc = {}
            desc['name'] = ''
            desc['low_limit'] = ''
            desc['producer'] = ''
            desc['supplyer'] = ''
            desc['price'] = ''
            desc['price_units'] = ''
            desc['price_date'] = datetime.datetime.now().date().strftime('%d.%m.%Y')
            desc['price_date_format'] = '%d.%m.%Y'
            desc['common_description'] = ''
            desc['mol_inchi'] = ''
            desc['mol_props'] = ''
            desc['mol_lumiprobe_data'] = ''
            desc['adduct_inchi'] = ''
            desc['adduct_props'] = ''
            self.data_obj['smart'] = desc
        else:
            try:
                self.data_obj = json.loads(data[4])
            except:
                self.data_obj = {
                'simple': {
                    'name': data[1],
                    'common_description': data[4],
                    'price_units': data[3],
                    'low_limit': data[5],
                           },
                'smart' : {}
                        }

        self.new = new
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                self.phyChemForm = phys_chem_props_interface()
                if len(self.data) > 1:
                    self.phyChemForm.unicode.value = self.data[2]
                if self.data_obj['smart'] == {}:
                    self.phyChemForm.description.value = self.data_obj['simple']['common_description']
                    self.phyChemForm.compound_name.value = self.data_obj['simple']['name']
                    self.phyChemForm.units.value = self.data_obj['simple']['price_units']
                    self.phyChemForm.compound_limit.value = self.data_obj['simple']['low_limit']
                else:
                    self.phyChemForm.description.value = self.data_obj['smart']['common_description']
                    self.phyChemForm.compound_name.value = self.data_obj['smart']['name']
                    self.phyChemForm.units.value = self.data_obj['smart']['price_units']
                    self.phyChemForm.compound_limit.value = self.data_obj['smart']['low_limit']
                    self.phyChemForm.producer.value = self.data_obj['smart']['producer']
                    self.phyChemForm.supplyer.value = self.data_obj['smart']['supplyer']
                    self.phyChemForm.price.value = self.data_obj['smart']['price']
                    self.phyChemForm.units.value = self.data_obj['smart']['price_units']
                    self.phyChemForm.structure.value = self.data_obj['smart']['mol_inchi']
                    if self.data_obj['smart']['mol_inchi'] != '':
                        self.phyChemForm.draw_structure()
                    if self.data_obj['smart']['mol_lumiprobe_data'] != '':
                        self.mol_desc = json.loads(self.data_obj['smart']['mol_lumiprobe_data'])
                        self.phyChemForm.mol_description.value = str(self.mol_desc)
                        self.phyChemForm.mol_lumiprobe_dict = json.loads(self.data_obj['smart']['mol_lumiprobe_data'])
                    self.phyChemForm.structure_adduct.value = self.data_obj['smart']['adduct_inchi']
        self.phyChemForm.on_save = self.do_save
        self.phyChemForm.on_cencel = self.do_cencel

    def on_save_data(self, data):
        pass

    def do_cencel(self):
        self.dialog.close()

    def do_save(self, data):
        self.data_obj['smart'] = data
        self.on_save_data(self.data_obj)
        self.dialog.close()


class add_widgets_dialog():
    def __init__(self, data, wis):
        self.data = data
        self.wis_list = []
        if len(wis) > 0:
            self.wis = json.loads(wis[0][4])
            self.wis_list = [key for key in self.wis.keys()]
        df = pd.DataFrame([i for i in self.data.values()])
        dict = {
            'name': list(df[1]),
            'unicode': list(df[2]),
            'label': ['' for i in range(df.shape[0])]
        }
        colDefs = [
            {"field": "name"},
            {"field": "label", 'editable': True},
        ]
        rowData = pd.DataFrame(dict)
        with ui.dialog() as self.dialog:
            with ui.card():
                self.selected_group = ui.label(text='')
                self.wis_group = ui.select(label='Установите группу:', options=self.wis_list, with_input=True,
                                            on_change=self.on_select_wis_event).classes('w-[400px]').style(
                    'font-size: 20px;')
                self.ag_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': rowData.to_dict('records'),
                        'rowSelection': 'multiple',
                        "pagination": True,
                        # "enableRangeSelection": True,
                    }
                    ,
                    theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
                self.ag_grid.auto_size_columns = True
                self.tab_rowdata = rowData.to_dict('records')
                self.ag_grid.on("cellValueChanged", self.update_cell_data)

                with ui.row():
                    self.new_group_name = ui.input(label='Имя новой группы').classes(
                        'w-[200px]').style('font-size: 20px;')
                    self.new_group_btn = ui.button('добавить группу', color='green',
                                                   on_click=self.on_new_group_event)

                with ui.row():
                    ui.button('Добавить', on_click=self.on_add_widgets_event)
                    ui.button('Отмена', on_click=self.dialog.close)

    def on_addjas_new_widgets(self, data):
        print(data)

    def on_new_group_event(self):
        self.group = self.new_group_name.value
        self.selected_group.set_text(self.group)

    def on_add_widgets_event(self):
        if self.group not in self.wis_list:
            self.wis[self.group] = []
        for row in self.tab_rowdata:
            self.wis[self.group].append([row['unicode'], row['label']])

        self.on_addjas_new_widgets(self.wis)
        self.dialog.close()

    def update_cell_data(self, e):
        self.tab_rowdata[e.args["rowIndex"]] = e.args["data"]

    def on_select_wis_event(self):
        self.group = self.wis_group.value
        self.selected_group.set_text(self.group)


class delete_widget_dialog():
    def __init__(self, data):
        self.data = data

        df = pd.DataFrame(self.data['group_data'])
        rowData = pd.DataFrame({
            'name': df[1],
            'delete': [False for i in range(df.shape[0])],
            'unicode': df[0]
        })

        colDefs = [
            {"field": "name"},
            {"field": "delete", 'editable': True},
        ]
        with ui.dialog() as self.dialog:
            with ui.card():
                self.ag_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': rowData.to_dict('records'),
                        'rowSelection': 'multiple',
                        "pagination": True,
                        # "enableRangeSelection": True,
                    }
                    ,
                    theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
                self.ag_grid.auto_size_columns = True
                self.tab_rowdata = rowData.to_dict('records')
                self.ag_grid.on("cellValueChanged", self.update_cell_data)

                with ui.row():
                    ui.button('Удалить', on_click=self.do_delete)
                    ui.button('Отмена', on_click=self.dialog.close)

    def update_cell_data(self, e):
        self.tab_rowdata[e.args["rowIndex"]] = e.args["data"]

    def on_delete(self, data, meta):
        print(data)

    def do_delete(self):
        edited = []
        for row in self.tab_rowdata:
            if not row['delete']:
                edited.append([row['unicode'], row['name']])
        self.data['wi_list'][self.data['group']] = edited
        self.on_delete(self.data['wi_list'], self.data['meta'])
        self.dialog.close()




class writeOff_dialog(api_db_interface):
    def __init__(self, api_IP, db_port, pincode, data, write_off=True):
        super().__init__(api_IP, db_port)

        self.write_off = write_off
        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.pincode = pincode
        self.data = data

        if self.write_off:
            text = 'Списание'
            title = 'Списать материал:'
            self.tab_name = 'output_tab'
        else:
            text = 'Поступление'
            title = 'Добавить материал:'
            self.tab_name = 'input_tab'

        if self.data != {}:
            self.unicode = self.data[0][2]

            if len(data) == 1:
                with ui.dialog() as self.dialog:
                    with ui.card():
                        ui.label(title)
                        ui.label(f'{self.data[0][1]}')
                        self.amount = ui.input()
                        with ui.row():
                            ui.button(text, on_click=self.on_write_off_click)
                            ui.button('Отмена', on_click=self.dialog.close)
            elif len(data) > 1:
                with ui.dialog() as self.dialog:
                    with ui.card():
                        ui.label(title)
                        #ui.label(f'{self.data[0][1]}')

                        df = pd.DataFrame(self.data)
                        rowData = pd.DataFrame({
                            'name': df[1],
                            'amount': [0. for i in range(df.shape[0])],
                            'unicode': df[2]
                        })

                        colDefs = [
                            {"field": "name", 'editable': True},
                            {"field": "amount", 'editable': True},
                        ]

                        self.ag_grid = ui.aggrid(
                            {
                                'columnDefs': colDefs,
                                'rowData': rowData.to_dict('records'),
                                'rowSelection': 'multiple',
                                "pagination": True,
                                # "enableRangeSelection": True,
                            }
                            ,
                            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
                        self.ag_grid.auto_size_columns = True
                        self.tab_rowdata = rowData.to_dict('records')
                        self.ag_grid.on("cellValueChanged", self.update_cell_data)

                        with ui.row():
                            ui.button(text, on_click=self.on_write_off_click_list)
                            ui.button('Отмена', on_click=self.dialog.close)

    def update_cell_data(self, e):
        self.tab_rowdata[e.args["rowIndex"]] = e.args["data"]
        #print(self.tab_rowdata)

    def on_write_off_click_list(self):
        self.substruct_from_stock(self.tab_name, self.tab_rowdata)
        #print(self.tab_rowdata)
        self.dialog.close()
        ui.run_javascript('location.reload();')

    def on_write_off_click(self):
        self.substruct_from_stock(self.tab_name,
                                  [{'name': self.data[0][1], 'amount': self.amount.value, 'unicode': self.data[0][2]}])
        self.dialog.close()
        ui.run_javascript('location.reload();')

    def get_all_data_in_tab_key(self, tab_name, key, value):
        url = f'{self.api_db_url}/get_keys_data/{self.db_name}/{tab_name}/{key}/{value}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()

    def substruct_solution(self, unicode, amount, tab_name):
        ret = self.get_all_data_in_tab_key('total_tab', 'unicode', unicode)
        data = json.loads(ret[0][4])
        if 'smart' in list(data.keys()):
            if data['smart']['mol_lumiprobe_data'] not in ['', '{}']:
                d_dict = json.loads(data['smart']['mol_lumiprobe_data'])
                for key in d_dict.keys():
                    amount_key = float(d_dict[key]) * float(amount)
                    r_ret = self.get_all_data_in_tab_key('total_tab', 'unicode', key)
                    url = f"{self.api_db_url}/insert_data/{self.db_name}/{tab_name}"
                    r = requests.post(url,
                                      json=json.dumps(
                                          [
                                              r_ret[0][1], key, amount_key,
                                              datetime.datetime.now().date().strftime(self.strftime_format),
                                              datetime.datetime.now().time().strftime(self.time_format),
                                              self.get_user_id()
                                          ]
                                      )
                                      , headers=self.headers())


    def substruct_from_stock(self, tab_name, rowdata):
        for row in rowdata:
            print(row)
            if row['name'].find('_sol_') > -1 and self.write_off:
                self.substruct_solution(row['unicode'], row['amount'], tab_name)
            else:
                if float(row['amount']) > 0:
                    url = f"{self.api_db_url}/insert_data/{self.db_name}/{tab_name}"
                    r = requests.post(url,
                    json=json.dumps(
                            [
                            row['name'], row['unicode'], row['amount'],
                            datetime.datetime.now().date().strftime(self.strftime_format),
                            datetime.datetime.now().time().strftime(self.time_format),
                            #user_id
                            self.get_user_id()
                            ]
                        )
                    , headers=self.headers())

    def get_user_id(self):
        url = f"{self.api_db_url}/get_keys_data/{self.db_name}/users/pin/{self.pincode}"
        r = requests.get(url, headers=self.headers())
        return r.json()[0][2]


class raw_mat_base_widget(api_db_interface):
    def __init__(self, api_IP, db_port, unicode, pincode):
        super().__init__(api_IP, db_port)

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.unicode = unicode
        self.pincode = pincode

    #[874, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-04', '12:02:17', '1783121115']
    #[879, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-10', '10:35:14', '1783121115']
    #[880, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-10', '10:36:06', '1783121115']

    def get_info_from_base(self):
        self.info_data = self.get_unicode_data_in_tab()
        self.remain = self.get_remaining_stock(self.unicode)
        self.output_data = self.get_unicode_output_data(self.unicode)
        #omap = oligomaps_search(self.db_IP, self.db_port)
        #map_data = omap.get_oligomaps_data()
        #print(map_data)
        #if len(map_data) > 0:
        #    self.accord_data = pd.DataFrame(map_data)['accord data']
        #    print(self.accord_data)

    def get_remaining_stock(self, unicode):
        url = f"{self.api_db_url}/get_remaining_stock/{self.db_name}/{unicode}"
        input_ret = requests.get(url, headers=self.headers())
        try:
            return input_ret.json()
        except:
            return {'exist': 0.}

    def get_unicode_data_in_tab(self):
        url = f'{self.api_db_url}/get_keys_data/{self.db_name}/total_tab/unicode/{self.unicode}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()

    def get_unicode_output_data(self, unicode):
        #url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/output_tab'
        #ret = requests.get(url, headers=self.headers())
        #pd.DataFrame(ret.json()).to_csv('output_tab.csv', sep='\t')
        url = f'{self.api_db_url}/get_keys_data/{self.db_name}/output_tab/unicode/{unicode}'
        ret = requests.get(url, headers=self.headers())

        return ret.json()

    def group_remain_unicode_data(self):
        df = pd.DataFrame(self.output_data)
        if df.shape[0] > 0:
            df['x_day'] = pd.to_datetime(df[4], format=self.strftime_format)
            df_g = df.groupby(pd.Grouper(key='x_day', freq='ME')).sum()
            df_g.reset_index(inplace=True)
            return list(df_g['x_day']), list(df_g[3])
        else:
            return [], []


class event():
    def __init__(self):
        self.type = 'null'
        self.x = 0
        self.y = 0

class menuItem():
    def __init__(self, name, color, canvas, x=10, y=10, w=110, h=40):
        self.name = name
        self.color = color
        self.canvas = canvas.add_layer()

        self.pos_x = x
        self.pos_y = y
        self.width = w
        self.height = h

        self.draw(event())

    def check_coord_click(self, e):
        if (e.image_x >= self.pos_x) and (e.image_x <= self.width + self.pos_x):
            if (e.image_y >= self.pos_y) and (e.image_y <= self.height + self.pos_y):
                return True
            else:
                return False
        else:
            return False

    def draw(self, e):
        fillop = 0.2
        if e.type == 'mousedown':
            fillop = 0.6

        self.canvas.content = ""
        self.canvas.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + 25}" '
                                      f'fill="white" font-size="20">{self.name}</text>')
        self.canvas.content += (f'<rect x={self.pos_x} y={self.pos_y} '
                               f'width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="{self.color}" fill-opacity="{fillop}"'
                                   f'stroke="{self.color}" stroke-width="4"/>')

    def on_click(self, e):
        pass

    def do_mousedown(self, e):
        if self.check_coord_click(e):
            self.draw(e)

    def do_mouseup(self, e):
        if self.check_coord_click(e):
            self.draw(e)
            self.on_click(e)


class infoPanel_menu():
    def __init__(self, img):
        self.items = {}
        self.image = img
        self.items['show info'] = menuItem('Show info', 'teal', self.image)
        self.items['write-off'] = menuItem('Write-off', 'orange', self.image, x=150)
        self.items['write-in'] = menuItem('Write-in', 'teal', self.image, x=270)
        self.items['add-btn'] = menuItem('Add reagent', 'teal', self.image, x=410, w=130)
        self.items['edit-btn'] = menuItem('Edit reagent', 'orange', self.image, x=550, w=130)
        self.items['add-widget-btn'] = menuItem('Add widgets', 'teal', self.image, x=710, w=130)

    def do_mousedown(self, e):
        for item in self.items.values():
            item.do_mousedown(e)

    def do_mouseup(self, e):
        for item in self.items.values():
            item.do_mouseup(e)

# red (красный)
# orange (оранжевый)
# amber (янтарный)
# yellow (жёлтый)
# lime (лаймовый)
# green (зелёный)
# emerald (изумрудный)
# teal (сине-зелёный)
# cyan (голубой)
# sky (небесный)
# blue (синий)
# indigo (индиго)
# violet (фиолетовый)
# purple (пурпурный)
# fuchsia (фуксия)
# pink (розовый)
# rose (роза)
# slate (сланцевый)
# gray (серый)
# zinc (цинковый)
# neutral (нейтральный)
# stone (каменный)

#Каждый цвет имеет оттенки с номерами 50, 100, 200… до 950, где 50 — самый светлый, 950 — самый тёмный. Например, для синего это могут быть классы text-blue-500 или bg-blue-200.


class rawMatList(api_db_interface):
    def __init__(self, api_IP, db_port, pincode, image, label_obj):
        super().__init__(api_IP, db_port)
        self.image = image
        self.label_obj = label_obj
        self.pincode = pincode
        self.visible = True

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.width = 1000
        self.height = 800
        self.search_top = 10
        self.scrol_width = 20
        self.menu_height = 100

        self.search_color = 'orange'
        self.selection_color = 'yellow'
        self.border_color = 'orange'
        self.selected_list_color = 'green'

        self.search_layer = self.image.add_layer()
        self.search_selection_layer = self.image.add_layer()
        self.selected_list_layer = self.image.add_layer()

        self.scroll_layer = self.image.add_layer()
        self.scroll_len = 30
        self.scroll_pos = 2
        self.scroll_yy = 2
        self.scroll_down = False
        self.scroll_border = self.height - self.menu_height - self.scroll_len - 5

        self.pushed_obj = {}
        self.selected_list = {}
        self.selected_row = {}

        self.main_tab_data_df = pd.DataFrame(self.get_all_data_in_tab('total_tab'))
        self.search_dataframe = self.main_tab_data_df
        self.search_list_index = 0

        self.on_mouse_down = self.on_down_event
        self.on_mouse_click = self.on_click_event
        self.draw_scroll()

    def set_visible(self, value):
        self.visible = value
        if not value:
            d = {}
            d['self.search_layer'] = self.search_layer.content
            d['self.search_selection_layer'] = self.search_selection_layer.content
            d['self.selected_list_layer'] = self.selected_list_layer.content
            d['self.scroll_layer'] = self.scroll_layer.content

            self.search_layer.content = ""
            self.search_selection_layer.content = ""
            self.selected_list_layer.content = ""
            self.scroll_layer.content = ""

            app.storage.user['raw_list_dict'] = d
        else:
            if 'raw_list_dict' in list(app.storage.user.keys()):
                d = app.storage.user.get('raw_list_dict')
                self.search_layer.content = d['self.search_layer']
                self.search_selection_layer.content = d['self.search_selection_layer']
                self.selected_list_layer.content = d['self.selected_list_layer']
                self.scroll_layer.content = d['self.scroll_layer']



    def do_mousedown(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            if self.check_scroll_catching(x, y):
                self.scroll_down = True
                self.scroll_yy = y
            self.on_mouse_down(e)

    def do_mousemove(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            if self.scroll_down:
                self.move_scroll(y)
            if self.search_layer.content != "":
                self.draw_search_selection(x, y)


    def do_mouseup(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            self.scroll_down = False

    def do_mouseclick(self, e):
        if self.visible:
            self.on_mouse_click(e)


    def check_seldown_event(self, key, y):
        for key in self.search_dict.keys():
            if (y >= key[1]) and (y <= key[3]):
                return True
        return False

    def check_selection_list_area(self, e):
        if e.image_y >= self.menu_height + 10:
            return True
        else:
            return False

    def on_down_event(self, e):
        if self.check_selection_list_area(e):
            if not e.ctrl:
                self.selected_list = {}
            if 'key' in list(self.selected_row.keys()):
                if self.check_seldown_event(self.selected_row['key'], e.image_y):
                    self.pushed_obj = self.selected_row.copy()
            if 1 in list(self.selected_row.keys()):
                self.label_obj.set_text(f"Info panel:  <<  {self.pushed_obj[1]}  >>")
                ui.notify(self.selected_row[1])

        self.draw_selected_list()

    def on_click_event(self, e):
        if 'key' in list(self.pushed_obj.keys()):
            self.selected_list[self.pushed_obj['key']] = self.pushed_obj.copy()
        if e.shift:
            pass
        self.draw_selected_list()

    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')
        self.draw_scroll()

    def draw_scroll(self):
        self.scroll_layer.content = ""

        self.scroll_layer.content = (f'<rect x={self.width - self.scrol_width} y={self.menu_height + self.scroll_pos} '
                                     f'width={15} height={self.scroll_len} '
                                     f' rx=10 ry=10 fill="none" '
                                     f'stroke="{self.border_color}" stroke-width="4"/>')

    def draw_search_layer(self, df):
        self.search_layer.content = ""
        dd = df.to_dict('records')
        x, y = 20, self.menu_height
        self.search_dict = {}
        for row in dd[self.search_list_index:]:
            d = row.copy()
            d['key'] = (x, y, self.width - 20, y + 35)
            self.search_dict[(x, y, self.width - 20, y + 35)] = d

            self.search_layer.content += (f'<text x="{20}" y="{y + 23}" '
                                          f'fill="white" font-size="20">{row[1]}</text>')

            self.search_layer.content += (f'<rect x={10} y={y} '
                                          f'width={self.width - 20 - self.scrol_width} height={35} '
                                          f'fill="{self.search_color}" fill-opacity="0.2"'
                                          f'stroke="{self.search_color}" stroke-width="2"/>')

            # print(y, self.height)
            y += 45
            if y > self.height:
                break

    def draw_selected_list(self):
        self.selected_list_layer.content = ""
        for row in self.selected_list.values():
            # print(row)
            key = row['key']
            self.selected_list_layer.content += (f'<rect x={10} y={key[1]} '
                                                 f'width={self.width - 20 - self.scrol_width} height={35} '
                                                 f'fill="{self.selected_list_color}" fill-opacity="0.2"'
                                                 f'stroke="{self.selected_list_color}" stroke-width="2"/>')

    def draw_search_selection(self, x, y):
        self.search_selection_layer.content = ""
        for key in self.search_dict.keys():
            if (y >= key[1]) and (y <= key[3]):
                self.selected_row = self.search_dict[key]
                self.search_selection_layer.content += (f'<rect x={10} y={key[1]} '
                                                        f'width={self.width - 20 - self.scrol_width} height={35} '
                                                        f'fill="{self.selection_color}" fill-opacity="0.2"'
                                                        f'stroke="{self.selection_color}" stroke-width="2"/>')

    def check_scroll_catching(self, x, y):
        if (x >= self.width - self.scrol_width) and (x <= self.width - self.scrol_width + 15):
            if (y >= self.menu_height + self.scroll_pos) and (
                    y <= self.menu_height + self.scroll_pos + self.scroll_len):
                return True
            else:
                return False
        else:
            return False

    def move_scroll(self, y):
        if (self.scroll_pos > 0) and (self.scroll_pos < self.scroll_border):
            self.scroll_pos = self.scroll_pos + y - self.scroll_yy
            if self.scroll_pos <= 0:
                self.scroll_pos = 2
            if self.scroll_pos >= self.scroll_border:
                self.scroll_pos = self.scroll_border - 1
            self.scroll_yy = y

        if self.search_dataframe.shape[0] >= 16:
            index = (self.search_dataframe.shape[0] - 16) * (self.scroll_pos - 2) / self.scroll_border
            index = int(round(index, 0))
        else:
            index = 0
        self.search_list_index = index
        self.draw_scroll()
        self.draw_search_layer(self.search_dataframe)

    def search_by_name(self, e):
        if e.value != '':
            self.search_dataframe = self.main_tab_data_df[self.main_tab_data_df[1].str.contains(e.value,
                                                                                                case=False, na=False)]
        else:
            self.search_dataframe = self.main_tab_data_df

        self.draw_search_layer(self.search_dataframe)


    def get_all_data_in_tab(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/{tab_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()


class descriptionPanel(raw_mat_base_widget):
    def __init__(self, api_IP, db_port, unicode, pincode, canvas, pos_x=20, pos_y=110, width=960, height=660):
        super().__init__(api_IP, db_port, unicode, pincode)

        self.main_layer = canvas
        self.visible = True

        self.pos_x = pos_x
        self.pos_y = pos_y
        self.width= width
        self.height = height
        self.color = 'neutral'#'orange' # neutral
        self.chart_color = 'orange'
        self.bar_color = 'green'
        self.chart_height = 250

        self.get_info_from_base()
        self.x_data, self.y_remain = self.group_remain_unicode_data()
        #print(self.info_data)
        #print(self.remain)
        #[[1, 'Ацетонитрил. ХЧ', 'INIT_BASE_CODE_OLIGO_LAB_0000001', '1L flask',
        #  'Ацетонитрил. ХЧ; 80 ppm H2O; для синтеза (необходимо сушить над ситами 3 или 4 ангстрема) и ВЭЖХ', 20, 1]]
        #{'exist': 21.0}
        self.draw()

    def draw_remain_monthly_data(self):
        main_x, main_y, main_w, main_h = self.pos_x + 10, self.pos_y + 45, self.width - 20, self.chart_height

        number = 12

        max_y_list = max(self.y_remain[-number:])
        self.avп_cons = round(sum(self.y_remain[-number:]) / number, 2)
        self.summ_cons = round(sum(self.y_remain[-number:]), 2)
        max_y = 190

        y_list = [int(round(i * max_y / max_y_list, 0)) for i in self.y_remain[-number:]]

        coord_x, x_width = 40, self.width // (number * 2)
        for x, y, data in zip(self.x_data[-number:], y_list, self.y_remain[-number:]):
            month_text = x.strftime(format="%m")
            self.main_layer.content += (f'<rect x={coord_x} y={main_h - y + self.pos_y + 20} '
                                        f'width={x_width} height={y} '
                                        f' rx=10 ry=10 fill="{self.bar_color}" fill-opacity="{0.6}"'
                                        f'stroke="{self.bar_color}" stroke-width="4"/>')
            self.main_layer.content += (f'<text x="{coord_x}" y="{main_h + self.pos_y + 40}" '
                                        f'fill="white" font-size="20">{month_text} m</text>')
            self.main_layer.content += (f'<text x="{coord_x}" y="{main_h - y + self.pos_y + 10}" '
                                        f'fill="white" font-size="20">{round(data, 2)}</text>')
            coord_x += x_width * 2

        self.main_layer.content += (f'<text x="{40}" y="{main_y + 25}" '
                                    f'fill="white" font-size="20"> средний расход: {self.avп_cons}</text>')
        self.main_layer.content += (f'<text x="{40}" y="{main_y + 45}" '
                                    f'fill="white" font-size="20"> расход за год: {self.summ_cons}</text>')

        self.main_layer.content += (f'<rect x={main_x} y={main_y} '
                                    f'width={main_w} height={main_h} '
                                    f' rx=10 ry=10 fill="{self.chart_color}" fill-opacity="{0.1}"'
                                    f'stroke="{self.chart_color}" stroke-width="4"/>')

    def draw_desc_string(self, text_line, pos_y):
        self.main_layer.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + pos_y}" '
                                f'fill="white" font-size="20">{text_line}</text>')


    def extract_text_desc(self):
        init_s = self.info_data[0][4]
        #print(init_s)
        try:
            dict = json.loads(init_s)
            if dict['smart'] != {}:
                if dict['smart']['mol_lumiprobe_data'] != '{}':
                    return str(json.loads(dict['smart']['mol_lumiprobe_data']))
                else:
                    return dict['smart']['common_description']
            else:
                return init_s
        except:
            return init_s

    def draw_description(self):
        #init_s = self.info_data[0][4]
        init_s = self.extract_text_desc()
        curs, str_curs = 0, 0
        text = ''
        pos_y = self.chart_height + 100
        drawed = False
        while True:
            space = init_s.find(' ', curs)
            if space == -1:
                space = len(init_s)
            if len(init_s[str_curs:space]) <= 70:
                text += f' {init_s[curs: space]}'
                drawed = False
            else:
                text += f' {init_s[curs: space]}'
                self.draw_desc_string(text, pos_y)
                drawed = True
                str_curs = space + 1
                text = ''
                pos_y += 25
            curs = space + 1
            if curs >= len(init_s):
                if not drawed and text != '':
                    self.draw_desc_string(text, pos_y)
                break

    def draw(self):
        fillop = 0.2

        self.main_layer.content = ""
        self.draw_description()
        self.draw_remain_monthly_data()

        remain = self.remain['exist']
        units = self.info_data[0][3]
        self.main_layer.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + 25}" '
                                f'fill="white" font-size="20">Остаток на складе: {round(remain, 0)} {units}</text>')
        self.main_layer.content += (f'<rect x={self.pos_x} y={self.pos_y} '
                                f'width={self.width} height={self.height} '
                                f' rx=10 ry=10 fill="{self.color}" fill-opacity="{fillop}"'
                                f'stroke="{self.color}" stroke-width="4"/>')

    def set_visible(self, value):
        self.visible = value
        if value:
            self.draw()
        else:
            self.main_layer.content = ""



class infoPanel(api_db_interface):

    def __init__(self, api_IP, db_port, pincode, lbl_obj):
        super().__init__(api_IP, db_port)
        self.api_IP, self.db_port = api_IP, db_port
        self.pincode = pincode
        self.label_obj = lbl_obj
        self._bkg = 'infopanel_bkg_2.png'

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"
        self.widget_db = 'gui_object_content_1.db'

        self.width = 1000
        self.height = 800
        self.search_top = 10
        self.scrol_width = 20
        self.menu_height = 100

        self.search_color = 'orange'
        self.selection_color = 'yellow'
        self.border_color = 'orange'
        self.selected_list_color = 'green'

        self.image = ui.interactive_image(f'/img/{self._bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])

        self.base_layer = self.image.add_layer()
        self.info_layer = self.image.add_layer()

        self.draw()
        self.info_menu = infoPanel_menu(self.image)
        self.rawMat = rawMatList(api_IP, db_port, pincode, self.image, self.label_obj)

        self.info_menu.items['show info'].on_click = self.on_show_info_menu_click
        self.info_menu.items['write-off'].on_click = self.on_write_off_stock
        self.info_menu.items['write-in'].on_click = self.on_write_in_stock
        self.info_menu.items['edit-btn'].on_click = self.on_edit_reagent_event
        self.info_menu.items['add-btn'].on_click = self.on_adjast_reagent_event
        self.info_menu.items['add-widget-btn'].on_click = self.on_add_widgets_event


    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')

    def on_show_info_menu_click(self, e):
        self.rawMat.set_visible(not self.rawMat.visible)

        if self.rawMat.visible:
            self.rawMat.draw_scroll()

        if not self.rawMat.visible:
            self.info_menu.items['show info'].color = 'red'
            self.info_menu.items['show info'].name = 'close info'

            if 2 in list(self.rawMat.pushed_obj.keys()):
                unicode = self.rawMat.pushed_obj[2]
            else:
                unicode = ''
            self.descript_panel = descriptionPanel(self.api_IP, self.db_port, unicode, self.pincode, self.info_layer)
        else:
            if self.descript_panel != None:
                self.descript_panel.set_visible(False)
            self.info_menu.items['show info'].color = 'teal'
            self.info_menu.items['show info'].name = 'show info'

    def on_show_info_menu_widget_click(self, data):

        self.rawMat.pushed_obj = data[0]
        self.rawMat.label_obj.set_text(f"Info panel:  <<  {self.rawMat.pushed_obj[1]}  >>")
        ui.notify(self.rawMat.pushed_obj[1])

        self.rawMat.set_visible(False)

        if not self.rawMat.visible:
            self.info_menu.items['show info'].color = 'red'
            self.info_menu.items['show info'].name = 'close info'

            unicode = data[0][2]

            self.descript_panel = descriptionPanel(self.api_IP, self.db_port, unicode, self.pincode, self.info_layer)
        else:
            if self.descript_panel != None:
                self.descript_panel.set_visible(False)
            self.info_menu.items['show info'].color = 'teal'
            self.info_menu.items['show info'].name = 'show info'


    def on_write_off_stock(self, e):
        if len(self.rawMat.selected_list.keys()) > 1:
            data = [i for i in self.rawMat.selected_list.values()]
            woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, data)
            woff.dialog.open()
            self.rawMat.selected_list = {}
        else:
            if self.rawMat.pushed_obj != {}:
                woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, [self.rawMat.pushed_obj])
                woff.dialog.open()
                self.rawMat.pushed_obj = {}
            else:
                ui.notify('Выберете объект со склада')


    def on_write_in_stock(self, e):
        if len(self.rawMat.selected_list.keys()) > 1:
            data = [i for i in self.rawMat.selected_list.values()]
            woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, data, write_off=False)
            woff.dialog.open()
            self.rawMat.selected_list = {}
        else:
            if self.rawMat.pushed_obj != {}:
                woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, [self.rawMat.pushed_obj],
                                       write_off=False)
                woff.dialog.open()
                self.rawMat.pushed_obj = {}
            else:
                ui.notify('Выберете объект со склада')


    def on_delete_wi_in_stock_event(self, data, meta):
        wi_list = {}
        for key in data.keys():
            if key != 'widgets':
                wi_list[key] = data[key]
        #print(meta)
        meta_list = [meta['#'], meta['object_id'], meta['date'], meta['date_format']]
        self.update_widgets([meta_list], wi_list)
        ui.run_javascript('location.reload();')

    def on_delete_wi_in_stock(self, data):
        del_dialog = delete_widget_dialog(data)
        del_dialog.on_delete = self.on_delete_wi_in_stock_event
        del_dialog.dialog.open()


    def update_tab(self, rowData):
        for row in rowData:
            url = f"{self.api_db_url}/update_data/{self.db_name}/total_tab/{row['#']}"
            ret = requests.put(url, json=json.dumps(
                {
                    'name_list': ['pos_name', 'unicode', 'units', 'description', 'lower_limit', 'actual'],
                    'value_list': [
                                            row['Name'], row['Unicode'],
                                            row['units'], row['Description'],
                                            row['low limit'], row['actual']
                    ]
                }
            ), headers=self.headers())


    def add_row(self, row):
        url = f"{self.api_db_url}/insert_data/{self.db_name}/total_tab"
        r = requests.post(url,
                          json=json.dumps(
                              [
                                  row['Name'], row['Unicode'], row['units'], row['Description'], row['low limit'], True
                              ]
                          )
                          , headers=self.headers())


    def on_edit_reagent_data(self, data):
        new_data = {}
        new_data['#'] = self.rawMat.pushed_obj[0]
        new_data['Name'] = data['smart']['name']
        new_data['Unicode'] = self.rawMat.pushed_obj[2]
        new_data['units'] = data['smart']['price_units']
        new_data['low limit'] = data['smart']['low_limit']
        new_data['Description'] = json.dumps(data)
        new_data['actual'] = self.rawMat.pushed_obj[6]
        #print(new_data)
        self.update_tab([new_data])
        ui.run_javascript('location.reload();')

    def on_edit_reagent_event(self, e):
        if self.rawMat.pushed_obj != {}:
            #print(self.rawMat.pushed_obj)
            ReagDialog = reagent_form_dialog(self.rawMat.pushed_obj)
            ReagDialog.on_save_data = self.on_edit_reagent_data
            ReagDialog.dialog.open()

    def gen_unicode(self):
        now = datetime.datetime.now()
        datetime_str = now.isoformat()
        hash_object = hashlib.sha256(datetime_str.encode('utf-8'))
        hash_hex = hash_object.hexdigest()
        return hash_hex

    def on_adjast_reagent_data(self, data):
        new_data = {}
        new_data['Name'] = data['smart']['name']
        new_data['Unicode'] = self.gen_unicode()
        new_data['units'] = data['smart']['price_units']
        new_data['low limit'] = data['smart']['low_limit']
        new_data['Description'] = json.dumps(data)
        new_data['actual'] = True
        #print(new_data)
        self.add_row(new_data)
        ui.run_javascript('location.reload();')

    def on_adjast_reagent_event(self, e):
        ReagDialog = reagent_form_dialog(self.rawMat.pushed_obj, new=True)
        ReagDialog.on_save_data = self.on_adjast_reagent_data
        ReagDialog.dialog.open()

    def load_widgets_content(self, object_id):
        url = f"{self.api_db_url}/get_keys_data/{self.widget_db}/main_tab/object_id/{object_id}"
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            return ret.json()
        else:
            return []

    def update_widgets(self, wis_meta, wi_list):
        row = {}
        row['#'] = wis_meta[0][0]
        row['object_id'] = wis_meta[0][1]
        row['date'] = wis_meta[0][2]
        row['date_format'] = wis_meta[0][3]
        row['obj_json'] = json.dumps(wi_list)

        #print(row)
        #print(self.headers())
        url = f"{self.api_db_url}/update_data/{self.widget_db}/main_tab/{row['#']}"
        ret = requests.put(url, json=json.dumps(
                {
                    'name_list': ['object_id', 'date', 'date_format', 'obj_json'],
                    'value_list': [
                                    row['object_id'], row['date'], row['date_format'], row['obj_json']
                    ]
                }
            ), headers=self.headers())

    def on_addjas_new_widgets(self, data):
        self.update_widgets(self.wis_data, data)
        ui.run_javascript('location.reload();')


    def on_add_widgets_event(self, e):
        self.wis_data = self.load_widgets_content('raw_material_widgets')
        add_wi_dialog = add_widgets_dialog(self.rawMat.selected_list, self.wis_data)
        add_wi_dialog.on_addjas_new_widgets = self.on_addjas_new_widgets
        add_wi_dialog.dialog.open()


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'mousedown':
            self.rawMat.do_mousedown(e)
            self.info_menu.do_mousedown(e)

        if e.type == 'mousemove':
            self.rawMat.do_mousemove(e)

        if e.type == 'mouseup':
            self.rawMat.do_mouseup(e)
            self.info_menu.do_mouseup(e)

        if e.type == 'click':
            self.rawMat.do_mouseclick(e)




class rawMatWidget(raw_mat_base_widget):

    def __init__(self, api_IP, db_port, unicode, pincode, on_click=lambda unicode: ui.notify(unicode), lbl='raw mat'):
        super().__init__(api_IP, db_port, unicode, pincode)

        self.label = lbl
        self.widget_bkg = 'widget_bkg_4.png'

        self.on_mouse_click = on_click

        self.width = 119
        self.height = 49
        self.title_height = 20

        self.title_color = 'orange'
        self.border_color = 'orange'

        self.image = ui.interactive_image(f'/img/{self.widget_bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])
        self.base_layer = self.image.add_layer()

        self.info_layer = self.image.add_layer()

        self.get_info_from_base()

        self.draw()
        self.draw_info()

    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=15 ry=15 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="2"/>')
        #self.base_layer.content += (f'<rect x=6 y=6 width={self.width - 12} height={self.title_height} rx=10 ry=10 '
        #                            f'fill="none" '
        #                           f'stroke="{self.title_color}" stroke-width="2"/>')

    def draw_info(self):
        self.info_layer.content = ""
        name = self.label
        self.info_layer.content += (f'<text x="{10}" y="{15}" '
                                                  f'fill="white" font-size="12">{name}</text>')
        fill_color = 'green'
        fill_const = 20
        fill_width_div = 1
        try:
            if self.remain['exist'] <= self.info_data[0][5]:
                fill_color = 'red'
                fill_width_div = 2
                fill_const = 0
        except:

            fill_color = 'red'
            fill_width_div = 2
            fill_const = 0
        #self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 25}" '
        #                            f'fill="white" font-size="10">{self.info_data[0][3]}</text>')
        amount = f"{round(self.remain['exist'], 2)} / {self.info_data[0][5]}"
        self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 15}" '
                                    f'fill="white" font-size="10">{amount}</text>')
        self.info_layer.content += (f'<rect x=8 y={self.title_height} '
                                    f'width={self.width // fill_width_div - fill_const} height={self.title_height} '
                                    f'fill="{fill_color}" fill-opacity="0.4"'
                                    f'stroke="{fill_color}" stroke-width="2"/>')


    def on_mouse_click(self, unicode):
        pass#print('click')


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'click':
            self.get_info_from_base()
            self.on_mouse_click(self.info_data)
