import pandas as pd
from nicegui import ui, app, events
from datetime import datetime
from datetime import timedelta
from OligoMap_utils import api_db_interface
import json
import requests
from oligoMass import molmassOligo as mmo
import asyncio
import pandas as pd
from molseq_lang import single_nucleic_acid_chain
from collections import Counter
import base64
from io import BytesIO


class price_db_manager(api_db_interface):
    def __init__(self):
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        super().__init__(IP, port)
        self.pincode = app.storage.user.get('pincode')
        self.db_price_name = 'oligo_price_blmx_2.db'

    def get_all_data_in_tab_key(self, tab_name, key, value):
        url = f'{self.api_db_url}/get_keys_data/{self.db_price_name}/{tab_name}/{key}/{value}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()

    def get_all_data_from_base(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_price_name}/{tab_name}'
        r = requests.get(url, headers=self.headers())
        if r.status_code == 200:
            return r.json()
        else:
            return []

    def get_all_types_groups(self):
        scale_list = self.get_all_scales()
        data = []
        for scale in scale_list:
            price_dict, type_dict = self.get_price_data(scale)
            for key in price_dict.keys():
                d = {}
                d['mod'] = key
                d['type'] = type_dict[key]
                data.append(d)
        df = pd.DataFrame(data)
        df = df.groupby('type').agg({'mod': list})
        df.reset_index(inplace=True)

        end5 = set(df[df['type']=='5end']['mod'].max())
        end3 = set(df[df['type']=='3end']['mod'].max())
        pt = set(df[df['type']=='purif_type']['mod'].max())
        return list(end5), list(end3), list(pt)


    def get_all_scales(self):
        tab = self.get_all_data_from_base('main')
        df = pd.DataFrame(tab)
        if df.shape[0] > 0:
            return list(df[1].unique())
        else:
            return []

    def get_price_data(self, scale):
        data = self.get_all_data_in_tab_key(tab_name='main', key='scale', value=scale)
        df = pd.DataFrame(data)
        if df.shape[0] > 0:
            df.sort_values(0, ascending=False, inplace=True)
            return json.loads(df.to_dict('records')[0][4]), json.loads(df.to_dict('records')[0][6])
        else:
            return {}

    def insert_price_data(self, scale, data_dict, type_dict):
        tab = self.get_all_data_from_base('main')
        insert = True
        for row in tab:
            if row[1] == scale:
                if row[4] == json.dumps(data_dict) and row[6] == json.dumps(type_dict):
                    insert = False
                else:
                    insert = True
        if insert:
            url = f'{self.api_db_url}/insert_data/{self.db_price_name}/main'
            r = requests.post(url,
                          json=json.dumps([
                              scale,
                              datetime.now().date().strftime('%d.%m.%Y'),
                              '%d.%m.%Y',
                              json.dumps(data_dict),
                              self.get_user_id(),
                              json.dumps(type_dict)
                          ]), headers=self.headers())
            if r.status_code == 200:
                ui.notify(f'Цены масштаба {scale} добавлены в базу')
        else:
            ui.notify(f'Данные уже есть в базе')

    def get_user_id(self):
        url = f"{self.api_db_url}/get_keys_data/{self.db_users}/users/pass/{self.pincode}"
        r = requests.get(url, headers=self.headers())
        return r.json()[0][1]

class oligo_price_calculator():
    def __init__(self):
        self.price_base = price_db_manager()
        self.error = {'unknown mod': '', 'unknown scale': '', 'unknown purif': ''}
        self.price = 0

    def check_scale(self, scale):
        scale_list = self.price_base.get_all_scales()
        if scale in scale_list:
            return True
        else:
            return False

    def check_modification(self, mod_list, scale, purif):
        price_data, type_data = self.price_base.get_price_data(scale)
        count = Counter(mod_list)
        sum_price = 0.
        for mod in count.keys():
            err_ctrl = True
            if mod is None:
                mod = ''
            if mod in price_data.keys():
                err_ctrl = False
                sum_price += price_data[mod] * count[mod]
            if err_ctrl and mod != 'none':
                self.error['unknown mod'] = mod
        err_ctrl = True
        for key in price_data.keys():
            if purif is None:
                purif = ''
            if key in purif:
                err_ctrl = False
                sum_price += price_data[key]
                break
        if err_ctrl:
            self.error['unknown purif'] = purif
        return sum_price


    def check_oligo_params(self, row_dict):
        oligo = single_nucleic_acid_chain(row_dict['Sequence'])
        mod_5end = row_dict["5'-end"]
        mod_3end = row_dict["3'-end"]
        scale = row_dict["Amount_OE"]
        purif = row_dict["Purification"]
        if self.check_scale(scale):
            mod_list = oligo.chain
            mod_list.extend([mod_5end, mod_3end])
            self.price = self.check_modification(mod_list, scale, purif)
        else:
            self.error['unknown scale'] = scale




class price_dialog():
    def __init__(self, scale):
        self.price_base = price_db_manager()
        self.price_data, self.type_data = self.price_base.get_price_data(scale)
        self.scale = scale
        self.rowdata = []
        self.init_()

    def init_(self):

        colDefs = [
            {"field": "Symbol", 'editable': True},
            {"field": "Price", 'editable': True},
            {"field": "Type", 'editable': True},
        ]

        scale_list = self.price_base.get_all_scales()

        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                with ui.row():
                    self.scale_param = ui.input(label='Set scale').style('width: 200px; font-size: 20px')
                    self.scale_selected = ui.select(options=scale_list, with_input=True, on_change=self.on_change_scale,
                                                    value='1-3').classes('w-[200px]')
                    ui.button('Add modification', on_click=self.on_add_modification)
                rowdata = []
                if self.price_data == {}:
                    self.scale_param.value = self.scale
                    rowdata = [
                        {'Symbol': 'A', 'Price': 35},
                        {'Symbol': 'C', 'Price': 35},
                        {'Symbol': 'G', 'Price': 35},
                        {'Symbol': 'T', 'Price': 35},
                    ]
                else:
                    self.scale_param.value = self.scale
                    for key in self.price_data.keys():
                        rowdata.append(
                            {
                            'Symbol': key,
                            'Price': self.price_data[key],
                            'Type': self.type_data[key],
                             }
                        )
                self.grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': rowdata,
                        'rowSelection': 'multiple',
                        "pagination": True,
                        "enterNavigatesVertically": True,
                        "enterNavigatesVerticallyAfterEdit": True,
                        "singleClickEdit": True,
                    },
                    theme='alpine-dark').style('width: 600px; height: 800px')
                self.grid.auto_size_columns = True
                self.grid.on("cellValueChanged", self.update_grid_cell_data)
                self.rowdata = self.grid.options['rowData']
                with ui.row():
                    ui.button('Сохранить', color='orange', on_click=self.on_save_data)
                    ui.button('Отмена', on_click=self.dialog.close)

    def on_save_data(self):
        price_dict, type_dict = {}, {}
        self.grid.options['rowData'] = self.rowdata
        self.grid.update()
        for row in self.grid.options['rowData']:
            price_dict[row['Symbol']] = row['Price']
            type_dict[row['Symbol']] = row['Type']
        self.price_base.insert_price_data(self.scale_param.value, price_dict, type_dict)
        self.dialog.close()

    def update_grid_cell_data(self, e):
        self.rowdata[e.args["rowIndex"]] = e.args["data"]

    def on_add_modification(self):
        self.rowdata.append({'Symbol': '', 'Price': 1})
        self.grid.options['rowData'] = self.rowdata
        self.grid.update()

    def on_change_scale(self, scale):
        self.price_data, self.type_data = self.price_base.get_price_data(scale.value)
        self.scale = scale.value
        self.scale_param.value = self.scale
        self.rowdata = []
        for key in self.price_data.keys():
            self.rowdata.append(
                {
                'Symbol': key,
                'Price': self.price_data[key],
                'Type': self.type_data[key],
                 }
            )
        self.grid.options['rowData'] = self.rowdata
        self.grid.update()

    

class price_tab():
    def __init__(self):
        self.synt_scale = '3-10'
        self.data, columnDefs = self.get_price_tab(self.synt_scale)
        self.tab = ui.aggrid({
                'columnDefs': columnDefs,
                'rowData': self.data.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
            }
            ,
            theme='alpine-dark').classes('h-[1000px]')
        self.rowdata = self.tab.options['rowData']
        self.tab.on("cellValueChanged", self.update_cell_data)

    def update_cell_data(self, e):
        self.rowdata[e.args["rowIndex"]] = e.args["data"]
        self.tab.options['rowData'] = self.rowdata
        self.tab.update()

    def on_change_scale(self, scale):
        self.synt_scale = scale.value
        self.data, columnDefs = self.get_price_tab(self.synt_scale)
        self.tab.options['rowData'] = self.data.to_dict('records')
        self.tab.update()

    def get_price_tab(self, scale):

        columnDefs = [
            {"field": 'Unit', "headerName": 'Unit', 'editable': True},
            {"field": 'Price, RUB', "headerName": 'Price, RUB', 'editable': True},
            {"field": 'Number', "headerName": 'Number', 'editable': True},
            {"field": 'Sum, RUB', "headerName": 'Sum, RUB', 'editable': True},
        ]

        unit_list = ['simple N', 'BHQ1', 'BHQ2', 'LNA', '6FAM', 'HEX', 'SIMA', 'R6G', 'ROX', 'TAMRA', 'Cy5', 'Cy3',
                     'HPLC']
        price_1_3_list = [30, 600, 600, 700, 1000, 1500, 1500, 1800, 1000, 1000, 1000, 1000, 1300]
        price_3_10_list = [40, 1000, 1000, 1200, 1800, 2000, 2000, 2200, 1800, 1800, 1800, 1800, 1300]
        price_10_15_list = [45, 1700, 1700, 2300, 2700, 3000, 3000, 3200, 2700, 2700, 2700, 2700, 1300]
        price_15_20_list = [50, 2400, 2400, 3300, 3700, 4000, 4000, 4200, 3700, 3700, 3700, 3700, 1300]
        price_30_40_list = [80, 4000, 4000, 4300, 4700, 5000, 5000, 5200, 4700, 4700, 4700, 4700, 1300]

        price_df = pd.DataFrame({'Unit': unit_list})
        if scale == '1-3':
            price_df['Price, RUB'] = price_1_3_list
        elif scale == '3-10':
            price_df['Price, RUB'] = price_3_10_list
        elif scale == '10-15':
            price_df['Price, RUB'] = price_10_15_list
        elif scale == '15-20':
            price_df['Price, RUB'] = price_15_20_list
        elif scale == '30-40':
            price_df['Price, RUB'] = price_30_40_list

        price_df['Number'] = [0 for i in range(price_df.shape[0])]
        price_df['Sum, RUB'] = [0 for i in range(price_df.shape[0])]

        return price_df, columnDefs




class input_order_page_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)
        self.init_pincode()
        self.db_name = 'scheduler_oligolab_2.db'

        self.price_base = price_db_manager()
        self.end5_list, self.end3_list, self.pt_list = self.price_base.get_all_types_groups()
        self.end5_list.append('none')
        self.end3_list.append('none')

        self.clip_js = '''
            (event) => {
                event.preventDefault();
                const text = event.clipboardData.getData('text');
                emitEvent('clipboard', text);
            }
        '''

        self.page_content()
        #self.add_context_menu()

        self.paste_area = ui.element('div').style('width: 1000px; height: 100px; border: 1px solid black;'
                                                  ).on('paste', js_handler=self.clip_js)



    def add_context_menu(self):
        with ui.context_menu() as self.context_menu:
            self.paste_clip = ui.menu_item('Вставить из буфера')

        self.paste_clip.on('click', js_handler=self.clip_js)




    def page_content(self):

        ui.on('clipboard', self.on_clipboard_event)

        df = pd.DataFrame(
            {
                "#": [],
                'Name': [],
                "5'-end": [],
                'Sequence': [],
                "3'-end": [],
                'Amount_OE': [],
                'Purification': [],
            }
        )

        columnDefs = [
            {"field": '#', "headerName": '#', 'editable': True},
            {"field": 'Name', "headerName": 'Name', 'editable': True},
            {"field": "5'-end", "headerName": '5-end', 'editable': True},
            {"field": 'Sequence', "headerName": 'Sequence', 'editable': True},
            {"field": "3'-end", "headerName": '3-end', 'editable': True},
            {"field": 'Amount_OE', "headerName": 'Amount', 'editable': True},
            {"field": 'Purification', "headerName": 'Purif type', 'editable': True},
            {"field": 'Price', "headerName": 'Price', 'editable': True},
            {"field": 'Error', "headerName": 'Error', 'editable': True},
        ]

        with ui.row():
            self.invoce = ui.input(label='Invoce name:').classes('w-[400px]')
            self.client = ui.input(label='Client:').classes('w-[400px]')
            self.price = ui.input(label='Price:').classes('w-[200px]')
            self.invoce_price_btn = ui.button('Add invoce', on_click=self.on_add_to_base
                                                  ).classes('w-[200px]')
            ui.button('Edit price base', color='orange', on_click=self.on_edit_price)
            with ui.column():
                self.invoce_price_btn = ui.button('Get price', on_click=self.culc_modif_count).classes('w-[200px]')
                self.culc_price = ui.input(label='Summary price').classes('w-[200px]')
            ui.upload(label='Загрузить форму',
                      on_upload=self.handle_upload).props().classes("max-w-full")#"accept=.mzdata.xml"

        with ui.row():
            with ui.dropdown_button('5 end mod', auto_close=True).classes('w-[200px]') as self.on_set_5end_mod:
                for key in self.end5_list:
                    ui.item(key, on_click=lambda n=key: self.on_set_5end_mod_event(n))
            #self.on_print_invoce_passport.props['color'] = 'orange'

            with ui.dropdown_button('3 end mod', auto_close=True).classes('w-[200px]') as self.on_set_3end_mod:
                for key in self.end3_list:
                    ui.item(key, on_click=lambda n=key: self.on_set_3end_mod_event(n))
            #self.on_print_invoce_passport.props['color'] = 'orange'

            with ui.dropdown_button('purif type', auto_close=True).classes('w-[200px]') as self.on_set_pt_mod:
                for key in self.pt_list:
                    ui.item(key, on_click=lambda n=key: self.on_set_pt_mod_event(n))
            #self.on_print_invoce_passport.props['color'] = 'orange'

        self.input_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
            }
            ,
            theme='alpine-dark').style('width: 2000px; height: 1000px')

            #self.price_grid = price_tab()

        #self.scale_selected.on_value_change(self.price_grid.on_change_scale)

        self.rowdata = self.input_tab.options['rowData']
        self.input_tab.on('paste', js_handler=self.clip_js)
        self.input_tab.auto_size_columns = True
        self.input_tab.on("cellValueChanged", self.update_cell_data)

    def update_cell_data(self, e):
        self.rowdata[e.args["rowIndex"]] = e.args["data"]
        self.input_tab.options['rowData'] = self.rowdata
        self.input_tab.update()

    def handle_upload(self, e: events.UploadEventArguments):
        data = e.content.read()
        self.b64_string = base64.b64encode(data).decode()
        print(self.b64_string)

        content = e.content.read()
        excel_io = BytesIO(content)
        df = pd.read_excel(excel_io, engine='openpyxl')
        print(df)

    def handle_clipboard(self, e):
        cname = {}
        cname['0'] = 'Name'
        cname['1'] = "5'-end"
        cname['2'] = 'Sequence'
        cname['3'] = "3'-end"
        cname['4'] = 'Amount_OE'
        cname['5'] = 'Purification'
        rowdata = []
        for row in e.args.split('\n'):
            insert_row = {}
            for i, col in enumerate(row.split('\t')):
                if col == '':
                    insert_row[cname[str(i)]] = 'none'
                else:
                    insert_row[cname[str(i)]] = col
            rowdata.append(insert_row)
        df = pd.DataFrame(rowdata)
        df.dropna(inplace=True)
        self.input_tab.options['rowData'] = df.to_dict('records')
        self.rowdata = self.input_tab.options['rowData']
        self.input_tab.update()

    def on_edit_price(self):
        dlg = price_dialog('1-3')
        dlg.dialog.open()

    async def on_set_5end_mod_event(self, mod):
        selrows = await self.input_tab.get_selected_rows()
        index_list = pd.DataFrame(selrows)['#'].to_list()
        for i, row in enumerate(self.input_tab.options['rowData']):
            if row['#'] in index_list:
                self.input_tab.options['rowData'][i]["5'-end"] = mod
        self.input_tab.update()

    async def on_set_3end_mod_event(self, mod):
        selrows = await self.input_tab.get_selected_rows()
        index_list = pd.DataFrame(selrows)['#'].to_list()
        for i, row in enumerate(self.input_tab.options['rowData']):
            if row['#'] in index_list:
                self.input_tab.options['rowData'][i]["3'-end"] = mod
        self.input_tab.update()

    async def on_set_pt_mod_event(self, mod):
        selrows = await self.input_tab.get_selected_rows()
        index_list = pd.DataFrame(selrows)['#'].to_list()
        for i, row in enumerate(self.input_tab.options['rowData']):
            if row['#'] in index_list:
                self.input_tab.options['rowData'][i]["Purification"] = mod
        self.input_tab.update()

    def on_clipboard_event(self, e):
        cname = {}
        cname['0'] = 'Name'
        cname['1'] = "5'-end"
        cname['2'] = 'Sequence'
        cname['3'] = "3'-end"
        cname['4'] = 'Amount_OE'
        cname['5'] = 'Purification'
        rowdata = []
        for row in e.args.split('\n'):
            insert_row = {}
            for i, col in enumerate(row.split('\t')):
                if col == '':
                    insert_row[cname[str(i)]] = 'none'
                else:
                    insert_row[cname[str(i)]] = col
            rowdata.append(insert_row)
        df = pd.DataFrame(rowdata)
        df.dropna(inplace=True)
        df['#'] = [i + 1 for i in range(df.shape[0])]

        self.input_tab.options['rowData'] = df.to_dict('records')
        self.rowdata = self.input_tab.options['rowData']
        self.input_tab.update()


    def on_paste_clipboard(self):
        pass


    def init_pincode(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        else:
            self.pincode = ''



    def process_data(self):
        self.total_price = 0


        self.df_tab = pd.DataFrame(self.input_tab.options['rowData'])
        self.price_tab = pd.DataFrame(self.price_grid.tab.options['rowData'])

        df = self.df_tab[self.df_tab['Sequence'] != '']

        len_list = []
        for seq in df['Sequence']:
            len_list.append(mmo.oligoNASequence(seq).size())
        df['Lenght'] = len_list

        self.df_tab.loc[self.df_tab['Sequence'] != '', 'Lenght'] = len_list

        sum_ = df['Lenght'].sum()
        price_ = self.price_tab[self.price_tab['Unit'] == 'simple N']['Price, RUB'].max()

        self.price_tab.loc[self.price_tab['Unit'] == 'simple N', 'Number'] = sum_
        self.price_tab.loc[self.price_tab['Unit'] == 'simple N', 'Sum, RUB'] = sum_ * float(price_)

        lbl_list = list(self.price_tab[self.price_tab['Unit'] != 'simple N']['Unit'])
        for lbl in lbl_list:
            sum_5 = df[df["5'-end"] == lbl].shape[0]
            price_5 = self.price_tab[self.price_tab['Unit'] == lbl]['Price, RUB'].max()

            sum_3 = df[df["3'-end"] == lbl].shape[0]
            price_3 = self.price_tab[self.price_tab['Unit'] == lbl]['Price, RUB'].max()

            sum_purif = df[df["Purification"] == lbl].shape[0]
            price_purif = self.price_tab[self.price_tab['Unit'] == lbl]['Price, RUB'].max()

            if sum_5 > 0:
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Number'] = sum_5
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Sum, RUB'] = sum_5 * float(price_5)

            if sum_3 > 0:
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Number'] = sum_3
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Sum, RUB'] = sum_3 * float(price_3)

            if sum_purif > 0:
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Number'] = sum_purif
                self.price_tab.loc[self.price_tab['Unit'] == lbl, 'Sum, RUB'] = sum_purif * float(price_purif)

        self.total_price = self.price_tab['Sum, RUB'].sum()
        self.culc_price.value = str(self.total_price)
        self.price_grid.tab.options['rowData'] = self.price_tab.to_dict('records')
        self.price_grid.tab.update()


    def culc_modif_count(self):
        rowdata = []
        for row in self.input_tab.options['rowData']:
            oligo = oligo_price_calculator()
            oligo.check_oligo_params(row)
            d = row.copy()
            d['Price'] = oligo.price
            out_err = {}
            for key in oligo.error.keys():
                if oligo.error[key] != '':
                    out_err[key] = oligo.error[key]
            d['Error'] = str(out_err)
            rowdata.append(d)
        df = pd.DataFrame(rowdata)
        self.culc_price.value = df['Price'].sum()
        self.input_tab.options['rowData'] = rowdata
        self.input_tab.update()



    def on_add_to_base(self):
        invoce = self.invoce.value
        client = self.client.value
        data = self.input_tab.options['rowData']
        price = self.price.value
        self.add_invoce_to_base(invoce, client, data, price)


    def add_invoce_to_base(self, invoce, client, data, price):
        url = f"{self.api_db_url}/insert_data/{self.db_name}/invoice_tab"
        params = json.dumps({'send': False, 'price': price})
        r = requests.post(url, json=json.dumps([invoce, client, params]), headers=self.headers())

        url = f"{self.api_db_url}/get_all_tab_data/{self.db_name}/invoice_tab"
        r = requests.get(url, headers=self.headers())
        id_list = [row[0] for row in r.json()]
        invoce_id = max(id_list)

        out_dataframe = pd.DataFrame(data)
        out_dataframe = out_dataframe[out_dataframe['Sequence'] != '']

        out_dataframe['order id'] = invoce_id
        out_dataframe['client id'] = invoce_id
        out_dataframe['input date'] = datetime.now().strftime('%m.%d.%Y')
        out_dataframe['output date'] = datetime.now().strftime('%m.%d.%Y')
        out_dataframe['status'] = 'in queue'
        out_dataframe['Sequence'] = out_dataframe['Sequence'].str.upper()
        out_dataframe.loc[out_dataframe["5'-end"] == '', "5'-end"] = 'none'
        out_dataframe.loc[out_dataframe["3'-end"] == '', "3'-end"] = 'none'

        len_list = []
        for seq in out_dataframe['Sequence']:
            len_list.append(mmo.oligoNASequence(seq).size())

        out_dataframe['lenght'] = len_list

        for (client_id, order_id, input_date, output_date, status,
             sequence, end5, end3, amount, purification, lenght, name) in zip(
            out_dataframe['client id'],
            out_dataframe['order id'],
            out_dataframe['input date'],
            out_dataframe['output date'],
            out_dataframe['status'],
            out_dataframe['Sequence'],
            out_dataframe["5'-end"],
            out_dataframe["3'-end"],
            out_dataframe['Amount_OE'],
            out_dataframe['Purification'],
            out_dataframe['lenght'],
            out_dataframe['Name'],
        ):
            url = f"{self.api_db_url}/insert_data/{self.db_name}/orders_tab"
            r = requests.post(url, json=json.dumps([client_id, order_id, input_date, output_date,
                               status, name, sequence, end5, end3, amount, purification, lenght]),
                              headers=self.headers())