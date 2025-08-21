import pandas as pd
from nicegui import ui, app
from io import BytesIO
from invoce_chart import invoceChart
from OligoMap_utils import api_db_interface
from OligoMap_utils import oligomaps_search
import requests
import json
from oligoMass import molmassOligo as mmo

from datetime import datetime
from datetime import timedelta


class navigation_menu(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        dark = ui.dark_mode(True)

        self.pincode = ''

        navigation = ui.row()
        with navigation:
            ui.link('Home', '/').style('font-size: 24px;')
            ui.link('Invoces', '/invoce_panel').style('font-size: 24px;')
            ui.link('Oligo synthesis', '/oligosynth_panel').style('font-size: 24px;')
            ui.link('Raw materials', '/rawmaterials_panel').style('font-size: 24px;')


        self.on_pincode_change = ui.input(label='pincode', placeholder='enter pincode',
                                          on_change=self.on_pincode_change_event)

        self.init_data()
        #ui.button('Show pincode', on_click=lambda: ui.notify(app.storage.user.get('pincode').value))

    def on_pincode_change_event(self, text):
        app.storage.user['pincode'] = text.value
        self.pincode = text.value
        if self.check_pincode().status_code == 200:
            ui.notify('pincode verified')

    def init_data(self):
        self.pincode = ''
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
            self.on_pincode_change.value = self.pincode
            if self.check_pincode().status_code == 200:
                ui.notify('pincode verified')
            else:
                ui.notify('Please, Enter pincode')
                self.on_pincode_change.value = ''


class invoice_page_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.pincode = ''
        self.db_name = 'scheduler_oligolab_2.db'
        self.status_hist_db = 'oligo_status_history_1.db'

        rowData = pd.DataFrame({
            '#': [0],
            'invoce': [''],
            'client': [''],
            'input date': [''],
            'out date': [''],
            'number': [0],
            'in queue%': [''],
            'synth%': [''],
            'purif%': [''],
            'formul%': [''],
            'fin%': [''],
            'product days': ['3'],
            'status': [''],
            'send': [True],
            'value P': ['']
        })

        colDefs = [
                        {"field": "#", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "invoce", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "client", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "input date", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "out date", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "number", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "in queue%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "synth%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "purif%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "formul%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "fin%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "product days", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "status", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "send", 'editable': True, 'filter': 'agSetColumnFilter'},
                        {"field": "value P", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                    ]

        summary_panel = ui.row()
        with summary_panel:
            col1 = ui.column()
            with col1:
                self.total_oligos_invoces = ui.circular_progress(min=0, max=1000, value=445,
                                                                 size='96px', color='orange')
                ui.label('Total oligos')
            col2 = ui.column()
            with col2:
                self.in_queue_oligos_invoces = ui.circular_progress(min=0, max=1000, value=100,
                                                                    size='96px', color='orange')
                ui.label('In queue')
            col3 = ui.column()
            with col3:
                self.synth_oligos_invoces = ui.circular_progress(min=0, max=1000, value=45,
                                                                 size='96px', color='orange')
                ui.label('Synthesis')
            col4 = ui.column()
            with col4:
                self.purif_oligos_invoces = ui.circular_progress(min=0, max=1000, value=65,
                                                                 size='96px', color='orange')
                ui.label('Purification')
            col5 = ui.column()
            with col5:
                self.formul_oligos_invoces = ui.circular_progress(min=0, max=1000, value=180,
                                                                  size='96px', color='orange')
                ui.label('Formulation')
            col6 = ui.column()
            with col6:
                self.fin_oligos_invoces = ui.circular_progress(min=0, max=1000, value=230,
                                                               size='96px', color='orange')
                ui.label('Finished')
            col7 = ui.column()
            with col7:
                self.total_price_oligos_invoces = ui.circular_progress(min=0, max=5000000, value=880000,
                                                                       size='96px', color='orange')
                ui.label('Price')

            with ui.column().classes('w-[1900px]'):
                self.invoceChart = invoceChart()

        self.ag_grid = ui.aggrid(
            {
                'columnDefs': colDefs,
                'rowData': rowData.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                #"enableRangeSelection": True,
            }
        ,
        theme='alpine-dark').classes('h-[800px]') # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True
        self.invoice_tab_rowdata = rowData.to_dict('records')
        self.ag_grid.on("cellValueChanged", self.update_invoce_cell_data)

        status_list = [
            'in queue',
            'synthesis',
            'purification',
            'formulation',
            'finished',
            'in progress',
            'wasted in progress',
            'not compleated in finished',
        ]

        with ui.row():
            self.on_load_button = ui.button('load invoces',
                                            on_click=self.on_load_button_event).classes('w-[200px]')

            self.on_show_actual_invoces_button = ui.button('Show actual',
                                            on_click=self.on_show_actual_invoces_button_event).classes('w-[200px]')

            self.on_show_invoces_content = ui.button('invoce content', on_click=self.on_show_invoces_content_event,
                                                     color='green').classes('w-[200px]')
            self.on_update_invoces_tab = ui.button('update tab',
                                                   on_click=self.on_update_invoces_tab_event).classes('w-[200px]')
            self.on_print_invoce_passport = ui.button('print passport',
                            on_click=self.on_print_invoce_passport_event, color='#FFA000').classes('w-[200px]')
            self.progressbar = ui.spinner(size='md', color='#FFA000')
            self.on_send_oligos_button = ui.button('Send selection to map', on_click=self.on_send_oligos_button_event,
                                                   color='green').classes('w-[200px]')
            self.on_show_oligos_synt_button = ui.button('Show synth number',
                                            on_click=self.on_show_oligos_synt_button_event).classes('w-[200px]')

        with ui.row().classes('gap-5'):
            with ui.column():
                self.on_show_by_status = ui.select(options=status_list, with_input=True,
                                                   on_change=self.on_show_by_status_event).classes('w-[200px]')
                self.input_date()
                self.on_print_orders_date_range = ui.button('Print data',
                                                on_click=self.on_print_orders_date_range_event).classes('w-[200px]')
            with ui.column().classes('w-[2400px]'):
                self.set_order_tab()

        self.progressbar.visible = False
        self.on_load_button.props('id="on_load_button"')

        self.init_data()


    def input_date(self):
        with ui.input('Start date') as self.start_date:
            with ui.menu().props('no-parent-event') as menu:
                with ui.date().bind_value(self.start_date):
                    with ui.row().classes('justify-end'):
                        ui.button('Close', on_click=menu.close).props('flat')
            with self.start_date.add_slot('append'):
                ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')

        with ui.input('End date') as self.end_date:
            with ui.menu().props('no-parent-event') as menu:
                with ui.date().bind_value(self.end_date):
                    with ui.row().classes('justify-end'):
                        ui.button('Close', on_click=menu.close).props('flat')
            with self.end_date.add_slot('append'):
                ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')

    def set_order_tab(self):

        invoce_content_df = pd.DataFrame(
            {
                '#': [1],
                'Name': [''],
                "5'-end": [''],
                'Sequence': [''],
                "3'-end": [''],
                'Amount, oe': ['5-10'],
                'Exist, oe': [0.],
                'Purification': ['Cart'],
                'Lenght': [''],
                'status': ['in queue'],
                'input date': [''],
                'output date': [''],
                'client id': [''],
                'order id': [''],
                'sufficiency': [0.],
                'synt, positions': [''],
            }
        )

        columnDefs = [
            {"field": "#", 'filter': 'agTextColumnFilter', 'floatingFilter': True,
             "checkboxSelection": True, "headerCheckboxSelection": True,},
            {"field": "Name", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "5'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Sequence", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "3'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Amount, oe", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Exist, oe", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Purification", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Lenght", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "status", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "input date", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "output date", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "order id", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "client id", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "synt, positions", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "sufficiency", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
        ]

        self.invoce_content_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': invoce_content_df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                ':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                                ' { background: "green" }',
            }
        ,
        theme='alpine-dark').classes('h-[1200px]') # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True

    def get_element_list(self, key):
        if key == 'button':
            return [
                ('on_load_button', 'click'),
                ('on_show_actual_invoces_button', 'click'),
                ('on_show_invoces_content', 'click'),
                ('on_print_invoce_passport', 'click'),
                ('on_send_oligos_button', 'click'),
                ('on_print_orders_date_range', 'click'),
                ('on_update_invoces_tab', 'click'),
                ('on_show_oligos_synt_button', 'click'),
            ]
        else:
            return []

    def save_passport(self, filename, data):
        buffer = BytesIO()
        data.to_excel(buffer, index=False)
        buffer.seek(0)
        ui.download(buffer.read(), filename=f'{filename}.xlsx')


    def update_invoce_cell_data(self, e):
        self.invoice_tab_rowdata[e.args["rowIndex"]] = e.args["data"]
        #print(e)
        self.on_update_invoces_tab.run_method('click')


    def get_all_invoces(self):
        url = f'{self.api_db_url}/get_all_invoces/{self.db_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()

    def get_invoce_content(self, selRows):
        out = pd.DataFrame(
            {
                '#': [1],
                'Name': [''],
                "5'-end": [''],
                'Sequence': [''],
                "3'-end": [''],
                'Amount, oe': ['5-10'],
                'Purification': ['Cart'],
                'Lenght': [''],
                'status': ['in queue'],
                'input date': [''],
                'output date': [''],
                'client id': [''],
                'order id': ['']
            }
        )
        out = out.to_dict('records')
        for row in selRows:
            out = []
            url = f"{self.api_db_url}/get_keys_data/{self.db_name}/orders_tab/order_id/{row['#']}"
            #orders_list = self.uny_db.get_all_tab_data_by_keys('orders_tab', 'order_id', row['#'])
            orders_list = requests.get(url, headers=self.headers())
            for r in orders_list.json():
                d = {}
                d['#'] = r[0]
                d['status'] = r[5]
                d['input date'] = r[3]
                d['output date'] = r[4]
                d['client id'] = row['client']
                d['order id'] = row['invoce']
                d["5'-end"] = str(r[8])
                d["Sequence"] = str(r[7])
                d["3'-end"] = str(r[9])
                d['Amount, oe'] = str(r[10])
                d['Purification'] = str(r[11])
                d['Lenght'] = str(r[12])
                d['Name'] = str(r[6])
                out.append(d)
            break
        return out

    def  get_oligo_prep_prognosis(self, d, cc = 1):
        data = self.get_invoce_content([d])
        df = pd.DataFrame(data)
        fam_count = df[df["5'-end"].str.contains('FAM')].shape[0]
        primer_count = df[df["5'-end"] == 'none'].shape[0]
        alk_count = df.shape[0] - fam_count - primer_count
        if primer_count == df.shape[0]:
            return 3 + cc
        elif primer_count <  df.shape[0] and alk_count == 0:
            return 5 + cc
        elif primer_count <  df.shape[0] and alk_count > 0:
            return 7 + cc

    def set_invoces_timing_prognosis(self, tab):
        out = []
        for row in tab:
            d = row.copy()
            n_days = self.get_oligo_prep_prognosis(d, cc=2)
            out_date = datetime.strptime(d['input date'], '%m.%d.%Y') + timedelta(days=n_days)

            start_date = datetime.strptime(d['input date'], '%m.%d.%Y')
            now_date = datetime.now()
            delta_days = (out_date - now_date).days
            prod_days = (now_date - start_date).days

            #if delta_days > 0:
            d['out date'] = out_date.strftime('%d.%m.%Y')

            if delta_days > 0:
                d['days_left'] = str(delta_days)
            elif delta_days == 0:
                d['days_left'] = '0+'
            else:
                d['days_left'] = '0'
            d['product days'] = f"{prod_days}/{n_days}"
            out.append(d)
        return out

    def get_summary_invoce_columns(self, sel_rowdata):
        out = {}
        if len(sel_rowdata) > 0:
            df = pd.DataFrame(sel_rowdata)
            key_list = ['number', 'in queue%', 'synth%', 'purif%', 'formul%', 'fin%', 'value P']
            for key in key_list:
                try:
                    df[key] = df[key].astype('int')
                    out[key] = str(df[key].sum())
                except:
                    out[key] = str(df[key].max())

                if key == 'fin%':
                    out[key] = sum([int(i[: i.find(' /')]) for i in list(df[key])])
        return out

    def on_show_actual_invoces_button_event(self):
        self.pincode = app.storage.user.get('pincode')

        df = pd.DataFrame(self.get_all_invoces())
        df = df[(df['status'] == 'in progress')|(df['send'] == False)]

        self.ag_grid.options['rowData'] = self.set_invoces_timing_prognosis(df.to_dict('records'))
        self.ag_grid.update()
        self.invoice_tab_rowdata = self.ag_grid.options['rowData']

        summary_dict = self.get_summary_invoce_columns(self.ag_grid.options['rowData'])
        self.total_oligos_invoces.value = int(summary_dict['number'])
        self.in_queue_oligos_invoces.value = int(summary_dict['in queue%'])
        self.synth_oligos_invoces.value = int(summary_dict['synth%'])
        self.purif_oligos_invoces.value = int(summary_dict['purif%'])
        self.formul_oligos_invoces.value = int(summary_dict['formul%'])
        self.fin_oligos_invoces.value = int(summary_dict['fin%'])
        self.total_price_oligos_invoces.value = int(summary_dict['value P'])

        app.storage.user['ag_grid_rowdata'] = self.ag_grid.options['rowData']
        app.storage.user['actual_invoces_summary_dict'] = summary_dict.copy()


    def on_load_button_event(self):
        self.pincode = app.storage.user.get('pincode')

        self.ag_grid.options['rowData'] = self.get_all_invoces()
        self.ag_grid.update()
        self.invoice_tab_rowdata = self.ag_grid.options['rowData']

        app.storage.user['ag_grid_rowdata'] = self.ag_grid.options['rowData']


    def load_history_ivoces(self):
        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        if len(self.pincode) == 4:
            url = f'{self.api_db_url}/get_all_tab_data/{self.status_hist_db}/main_tab'
            ret = requests.get(url, headers=self.headers())

            self.history_stat_data = []
            for row in ret.json():
                d = {}
                d['#'] = row[0]
                d['Date'] = row[1]
                d['d_format'] = '%d.%m.%Y'
                d['Time'] = row[2]
                d['data'] = json.loads(row[3])
                self.history_stat_data.append(d)
            self.invoceChart.draw_stat_data(self.history_stat_data)


    def update_send_invoce_data(self, rowdata):
        param_list = []
        for row in rowdata:
            param_list.append({'id': row['#'], 'send_param': json.dumps({'send': row['send'],
                                                                         'price': row['value P'],})})
        if len(param_list) > 0:
            url = f"{self.api_db_url}/send_invoces_update/{self.db_name}"
            #print(json.dumps(param_list))
            ret = requests.put(url, json=json.dumps(param_list), headers=self.headers())


    def on_update_invoces_tab_event(self):
        self.pincode = app.storage.user.get('pincode')

        self.update_send_invoce_data(self.invoice_tab_rowdata)
        df = pd.DataFrame(self.invoice_tab_rowdata)
        df = df[(df['status'] == 'in progress') | (df['send'] == False)]
        self.invoice_tab_rowdata = df.to_dict('records')
        self.ag_grid.options['rowData'] = self.invoice_tab_rowdata
        self.ag_grid.update()
        app.storage.user['ag_grid_rowdata'] = self.ag_grid.options['rowData']


    async def on_show_invoces_content_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.progressbar.visible = True
        selrows = await self.ag_grid.get_selected_rows()
        orders_tab = self.get_invoce_content(selrows)

        out = []
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        date = datetime.strptime(orders_tab[0]['input date'], '%m.%d.%Y')
        date = date - timedelta(days=10)
        oserch.map_list = pd.DataFrame(oserch.get_oligomaps_date_tail(date.strftime('%Y-%m-%d'), tail_len=30))
        for row in orders_tab:
            d = row.copy()
            maps = oserch.find_amount_by_order_id(row['#'])
            df = maps['Dens, oe/ml'] * maps['Vol, ml']
            d['Exist, oe'] = round(df.sum(), 0)
            limit = oserch.get_low_amount_limit(d['Amount, oe'])
            d['sufficiency'] = d['Exist, oe'] - limit
            out.append(d)

        self.invoce_content_tab.options['rowData'] = out
        self.invoce_content_tab.update()

        app.storage.user['invoce_content_tab_rowdata'] = self.invoce_content_tab.options['rowData']

        self.progressbar.visible = False


    def print_pass(self, rowData):
        out_tab = []
        index_ = 1

        for row in rowData:
            nseq = row['Seq']
            if type(row['Purif type']) == str:
                if row['Purif type'].find('_') > 0:
                    fluoro = row['Purif type'][row['Purif type'].find('_')+1:]
                    nseq = row['Seq'].replace('[Alk]', f"[{fluoro}]")
                    #print(nseq)
                    o = mmo.oligoNASequence(nseq)
                else:
                    o = mmo.oligoNASequence(row['Seq'])
            d = {}
            d['#'] = index_
            index_ += 1
            #d['Position'] = row['Position']
            d['Name'] = row['Name'] + f"  ({row['Synt number']}_{row['Position']})"
            d['Sequence'] = nseq
            d['Amount,_oe'] = int(round(row['Dens, oe/ml'] * row['Vol, ml'], 0))
            if o.getExtinction() > 0:
                d['Amount,_nmol'] = int(round(d['Amount,_oe'] * 1e6 / o.getExtinction(), 0))
            else:
                d['Amount,_nmol'] = 0.
            d['Desolving'] = int(d['Amount,_nmol'] * 10)

            d['Purification'] = row['Purif type']
            d['order_ID'] = row['Order id']
            d['Status'] = row['Status']
            try:
                d['Mass,_Da'] = round(o.getAvgMass(), 2)
            except:
                d['Mass,_Da'] = 'unknown modiff'
            d['Extinction'] = o.getExtinction()

            out_tab.append(d)
        return pd.DataFrame(out_tab)


    async def on_print_invoce_passport_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.progressbar.visible = True
        selrows = await self.ag_grid.get_selected_rows()
        orders_tab = self.get_invoce_content(selrows)
        print('PASSPORT')
        out = []
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        date = datetime.strptime(orders_tab[0]['input date'], '%m.%d.%Y')
        date = date - timedelta(days=10)
        oserch.map_list = pd.DataFrame(oserch.get_oligomaps_date_tail(date.strftime('%Y-%m-%d'), tail_len=30))
        for row in orders_tab:
            d = row.copy()
            maps = oserch.find_amount_by_order_id(row['#'])
            maps = maps[maps['Dens, oe/ml'] > 0]
            maps['Synt number'] = pd.to_numeric(maps['Synt number'], errors='coerce')
            maps = maps.sort_values(by='Synt number', ascending=False)
            #print(maps[['Synt number', 'Position', 'Dens, oe/ml', 'Vol, ml']])
            tab = maps.to_dict('records')
            if len(tab) > 0:
                df = maps['Dens, oe/ml'] * maps['Vol, ml']
                d['Exist, oe'] = round(df.sum(), 0)
                limit = oserch.get_low_amount_limit(d['Amount, oe'])
                d['sufficiency'] = d['Exist, oe'] - limit
                d['Purif type'] = maps['Purif type'].max()
                #d['Position'] = maps['Position'].min()
                d['Position'] = tab[0]['Position']
                #maps['Synt number'] = pd.to_numeric(maps['Synt number'], errors='coerce')
                #d['Synt number'] = maps['Synt number'].max()
                d['Synt number'] = tab[0]['Synt number']
                d['Status'] = maps['Status'].max()
                d['Vol, ml'] = 1
                d['Dens, oe/ml'] = d['Exist, oe']
                d['Order id'] = maps['Order id'].max()
                d['Seq'] = maps['Sequence'].max()
                out.append(d)

        self.invoce_content_tab.options['rowData'] = out
        self.invoce_content_tab.update()

        app.storage.user['invoce_content_tab_rowdata'] = self.invoce_content_tab.options['rowData']

        pass_filename = selrows[0]['invoce'].replace('/', '_')
        pass_data = self.print_pass(self.invoce_content_tab.options['rowData'])
        self.save_passport(pass_filename, pass_data)
        self.progressbar.visible = False


    async def on_show_oligos_synt_button_event(self):
        self.pincode = app.storage.user.get('pincode')
        selrows = await self.invoce_content_tab.get_selected_rows()
        sel_orders = list(pd.DataFrame(selrows)['#'])
        rowdata = self.invoce_content_tab.options['rowData']

        out = []
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        date = datetime.strptime(selrows[0]['input date'], '%m.%d.%Y')
        date = date - timedelta(days=10)
        oserch.map_list = pd.DataFrame(oserch.get_oligomaps_date_tail(date.strftime('%Y-%m-%d'), tail_len=30))
        for row in rowdata:
            d = row.copy()
            if row['#'] in sel_orders:
                maps = oserch.find_amount_by_order_id(row['#'])
                s = ""
                for syn, pos in zip(maps['Synt number'], maps['Position']):
                    #print(syn, pos)
                    s += f'{syn}_{pos}, '
                d['synt, positions'] = s
            out.append(d)

        self.invoce_content_tab.options['rowData'] = out
        self.invoce_content_tab.update()
        app.storage.user['invoce_content_tab_rowdata'] = self.invoce_content_tab.options['rowData']


    async def on_send_oligos_button_event(self):
        self.pincode = app.storage.user.get('pincode')
        selrows = await self.invoce_content_tab.get_selected_rows()

        if 'oligos_list_synth' in list(app.storage.user.keys()):
            app.storage.user['oligos_list_synth'].extend(selrows)
        else:
            app.storage.user['oligos_list_synth'] = []
            app.storage.user['oligos_list_synth'].extend(selrows)


    def get_orders_by_status(self, status):

        def get_in_progress(find_list = ['synthesis', 'purification', 'formulation']):
            out = []
            for st in find_list:
                url = f'{self.api_db_url}/get_orders_by_status/{self.db_name}/{st}'
                ret = requests.get(url, headers=self.headers())
                out.extend(ret.json())
            return out
        out = []
        if status == 'in progress':
            out = get_in_progress(find_list = ['synthesis', 'purification', 'formulation'])
        elif status == 'total data':
            out = get_in_progress(find_list = ['in queue', 'synthesis', 'purification', 'formulation', 'finished'])
        else:
            url = f'{self.api_db_url}/get_orders_by_status/{self.db_name}/{status}'
            ret = requests.get(url, headers=self.headers())
            out.extend(ret.json())
        return out


    def on_show_by_status_event(self, param):
        self.pincode = app.storage.user.get('pincode')
        out = self.get_orders_by_status(param.value)
        self.invoce_content_tab.options['rowData'] = out
        self.invoce_content_tab.update()
        app.storage.user['invoce_content_tab_rowdata'] = self.invoce_content_tab.options['rowData']


    def on_print_orders_date_range_event(self):
        self.pincode = app.storage.user.get('pincode')
        start_date = datetime.strptime(self.start_date.value, "%Y-%m-%d")
        end_date = datetime.strptime(self.end_date.value, "%Y-%m-%d")

        data = pd.DataFrame(self.invoce_content_tab.options['rowData'])
        data['Date'] = pd.to_datetime(data['input date'], format='%m.%d.%Y')
        conditions = (data['Date'] >= start_date)&(data['Date'] <= end_date)
        columns = ["#", "Name", "5'-end", "Sequence", "3'-end", "Amount, oe", "Purification", "Lenght", "order id",
                   "client id", "input date"]
        df = data[conditions][columns]
        self.save_passport('data_range', df)
        app.storage.user['invoce_content_tab_rowdata'] = self.invoce_content_tab.options['rowData']


    def init_data(self):
        self.ag_grid.options['rowData'] = app.storage.user.get('ag_grid_rowdata')
        self.ag_grid.update()
        self.load_history_ivoces()

        self.invoce_content_tab.options['rowData'] = app.storage.user.get('invoce_content_tab_rowdata')
        self.invoce_content_tab.update()

        try:
            summary_dict = app.storage.user.get('actual_invoces_summary_dict')
            self.total_oligos_invoces.value = int(summary_dict['number'])
            self.in_queue_oligos_invoces.value = int(summary_dict['in queue%'])
            self.synth_oligos_invoces.value = int(summary_dict['synth%'])
            self.purif_oligos_invoces.value = int(summary_dict['purif%'])
            self.formul_oligos_invoces.value = int(summary_dict['formul%'])
            self.fin_oligos_invoces.value = int(summary_dict['fin%'])
            self.total_price_oligos_invoces.value = int(summary_dict['value P'])
        except:
            pass

