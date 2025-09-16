import pandas as pd
from nicegui import ui, app
from datetime import datetime
from datetime import timedelta
from OligoMap_utils import api_db_interface
import json
import requests
from oligoMass import molmassOligo as mmo
import asyncio
import pandas as pd


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

        #self.paste_area = ui.element('div').style('width: 1000px; height: 100px; border: 1px solid black;'
        #                                          ).on('paste', js_handler='''
        #    (event) => {
        #        event.preventDefault();
        #        const text = event.clipboardData.getData('text');
        #        emitEvent('clipboard', text);
        #    }
        #''')


    def add_context_menu(self):
        with ui.context_menu() as self.context_menu:
            self.paste_clip = ui.menu_item('Вставить из буфера')

        self.paste_clip.on('click', js_handler=self.clip_js)




    def page_content(self):

        ui.on('clipboard', self.on_clipboard_event)

        df = pd.DataFrame(
            {
                'Name': [''],
                "5'-end": [''],
                'Sequence': [''],
                "3'-end": [''],
                'Amount_OE': ['5-10'],
                'Purification': ['Cart'],
            }
        )

        columnDefs = [
            {"field": 'Name', "headerName": 'Name', 'editable': True},
            {"field": "5'-end", "headerName": '5-end', 'editable': True},
            {"field": 'Sequence', "headerName": 'Sequence', 'editable': True},
            {"field": "3'-end", "headerName": '3-end', 'editable': True},
            {"field": 'Amount_OE', "headerName": 'Amount', 'editable': True},
            {"field": 'Purification', "headerName": 'Purif type', 'editable': True},
        ]

        scale_list = [
            '1-3',
            '3-10',
            '10-15',
            '15-20',
            '30-40',
        ]


        with ui.grid(columns=2).classes("w-full").style("grid-template-columns: 1500px 1000px"):

            with ui.row():
                self.invoce = ui.input(label='Invoce name:').classes('w-[400px]')
                self.client = ui.input(label='Client:').classes('w-[400px]')
                self.price = ui.input(label='Price:').classes('w-[200px]')
                self.invoce_price_btn = ui.button('Add invoce', on_click=self.on_add_to_base
                                                  ).classes('w-[200px]')

            with ui.row():
                self.scale_selected = ui.select(options=scale_list, with_input=True,
                value='3-10').classes('w-[400px]')
                self.invoce_price_btn = ui.button('Culc price', on_click=self.process_data).classes('w-[200px]')
                self.culc_price = ui.input(label='Culc price').classes('w-[200px]')

            self.input_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
            }
            ,
            theme='alpine-dark').classes('h-[1000px]')

            self.price_grid = price_tab()

        self.scale_selected.on_value_change(self.price_grid.on_change_scale)

        self.rowdata = self.input_tab.options['rowData']
        self.input_tab.on('paste', js_handler=self.clip_js)
        #self.input_tab.on('keydown.ctrl.v', js_handler='''
        #                        async () => emitEvent("clipboard", await navigator.clipboard.readText())
        #                        ''')
        self.input_tab.auto_size_columns = True
        self.input_tab.on("cellValueChanged", self.update_cell_data)

    def update_cell_data(self, e):
        self.rowdata[e.args["rowIndex"]] = e.args["data"]
        #print(self.rowdata[e.args["rowIndex"]].keys())
        #print(e.args["data"])
        self.input_tab.options['rowData'] = self.rowdata
        self.input_tab.update()
        #print(self.input_tab.options['rowData'][0].keys())

    def handle_clipboard(self, e):
        print(e)
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