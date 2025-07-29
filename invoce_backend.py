import json

import pandas as pd
import requests

from datetime import datetime
from datetime import timedelta

from OligoMap_utils import api_db_interface
from OligoMap_utils import oligomaps_search

from oligoMass import molmassOligo as mmo
import matplotlib.pyplot as plt


class invoce_table(api_db_interface):

    def __init__(self, api_IP, db_port, stack):
        super().__init__(api_IP, db_port)

        self.oligomap_stack = stack
        self.pincode = ''
        self.db_name = 'scheduler_oligolab_2.db'
        self.status_hist_db = 'oligo_status_history_1.db'
        self.client = {}
        self.client_frontend = {}


    def init_frontend(self, front, ip):
        if ip not in list(self.client.keys()):
            self.client[ip] = ''
            self.client_frontend[ip] = front.get_model()
        self.frontend = front
        self.frontend.set_model(self.client_frontend[ip])
        #self.frontend.on_pincode_change.value = self.client[ip]
        if self.client[ip] != '':
            self.load_history_ivoces()


    def init_application(self, app):
        self.app = app


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


    def get_all_invoces(self):
        url = f'{self.api_db_url}/get_all_invoces/{self.db_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()


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


    def get_pincode(self, ip):
        if ip in list(self.client.keys()):
            return self.client[ip]
        else:
            return ''


    def on_load_button(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.frontend.ag_grid.options['rowData'] = self.get_all_invoces()
        self.frontend.ag_grid.update()
        self.frontend.invoice_tab_rowdata = self.frontend.ag_grid.options['rowData']

        #df = pd.DataFrame(self.frontend.ag_grid.options['rowData'])
        #df['date group dt'] = pd.to_datetime(df['input date'], format='%m.%d.%Y')
        #df['date group'] = df['date group dt'].dt.strftime('%m-%Y')
        #df['invoce num'] = 1
        #df = df.groupby('date group').agg({
        #                                    'invoce num': 'sum',
        #                                    'number': 'sum',
        #                                    'date group dt': 'min'
        #                                   })
        #df.sort_values(by='date group dt', ascending=True, inplace=True)
        #df = df.reset_index()
        #print(df)
        #plt.style.use('dark_background')
        #plt.rcParams["figure.figsize"] = (18, 4)
        #plt.plot(df['date group'], df['number'], 'o-')
        #plt.xlabel('month')
        #plt.ylabel('Number of oligos')
        #plt.savefig('images/number_oligos_plot_1.png')

        #self.frontend.invoces_stat_mage.force_reload()

        self.client_frontend[ip] = self.frontend.get_model()


    def on_show_actual_invoces_button(self):
        print(self.client)
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]#self.frontend.input_pincode.value
        #print(self.pincode)

        df = pd.DataFrame(self.get_all_invoces())
        df = df[(df['status'] == 'in progress')|(df['send'] == False)]

        self.frontend.ag_grid.options['rowData'] = self.set_invoces_timing_prognosis(df.to_dict('records'))
        self.frontend.ag_grid.update()
        self.frontend.invoice_tab_rowdata = self.frontend.ag_grid.options['rowData']

        summary_dict = self.get_summary_invoce_columns(self.frontend.ag_grid.options['rowData'])
        #print(summary_dict)
        #print(int(summary_dict['number']))
        self.frontend.total_oligos_invoces.value = int(summary_dict['number'])
        self.frontend.in_queue_oligos_invoces.value = int(summary_dict['in queue%'])
        self.frontend.synth_oligos_invoces.value = int(summary_dict['synth%'])
        self.frontend.purif_oligos_invoces.value = int(summary_dict['purif%'])
        self.frontend.formul_oligos_invoces.value = int(summary_dict['formul%'])
        self.frontend.fin_oligos_invoces.value = int(summary_dict['fin%'])
        self.frontend.total_price_oligos_invoces.value = int(summary_dict['value P'])

        self.client_frontend[ip] = self.frontend.get_model()


    def update_send_invoce_data(self, rowdata):
        param_list = []
        for row in rowdata:
            param_list.append({'id': row['#'], 'send_param': json.dumps({'send': row['send'],
                                                                         'price': row['value P'],})})
        if len(param_list) > 0:
            url = f"{self.api_db_url}/send_invoces_update/{self.db_name}"
            #print(json.dumps(param_list))
            ret = requests.put(url, json=json.dumps(param_list), headers=self.headers())

    def on_update_invoces_tab(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.update_send_invoce_data(self.frontend.invoice_tab_rowdata)
        df = pd.DataFrame(self.frontend.invoice_tab_rowdata)
        df = df[(df['status'] == 'in progress') | (df['send'] == False)]
        self.frontend.invoice_tab_rowdata = df.to_dict('records')
        self.frontend.ag_grid.options['rowData'] = self.frontend.invoice_tab_rowdata
        self.frontend.ag_grid.update()

        self.client_frontend[ip] = self.frontend.get_model()


    async def on_show_invoces_content(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        self.frontend.progressbar.visible = True
        selrows = await self.frontend.ag_grid.get_selected_rows()
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

        self.frontend.invoce_content_tab.options['rowData'] = out
        self.frontend.invoce_content_tab.update()
        self.client_frontend[ip] = self.frontend.get_model()
        self.frontend.progressbar.visible = False


    def on_show_by_status(self, param):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        out = self.get_orders_by_status(param.value)
        #print(out[0].keys())

        self.frontend.invoce_content_tab.options['rowData'] = out
        #print(self.frontend.invoce_content_tab.options['rowData'][0].keys())
        self.frontend.invoce_content_tab.update()
        self.client_frontend[ip] = self.frontend.get_model()


    def on_pincode_change(self, text):
        self.client[self.app.storage.user.get('client_ip')] = text.value
        self.client_frontend[self.app.storage.user.get('client_ip')] = self.frontend.get_model()
        #self.frontend.ui.run_javascript('simulateClick()')
        self.load_history_ivoces()


    async def on_print_invoce_passport(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        self.frontend.progressbar.visible = True
        selrows = await self.frontend.ag_grid.get_selected_rows()
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
            maps = maps[maps['Dens, oe/ml'] > 0]
            maps['Synt number'] = pd.to_numeric(maps['Synt number'], errors='coerce')
            maps = maps.sort_values(by='Synt number', ascending=False)
            #print(maps[['Synt number', 'Position', 'Dens, oe/ml', 'Vol, ml']])
            tab = maps.to_dict('records')
            #print(tab[0])
            #print(maps[['Synt number', 'Position', 'Dens, oe/ml', 'Vol, ml']])
            #print(tab[0])
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

        self.frontend.invoce_content_tab.options['rowData'] = out
        self.frontend.invoce_content_tab.update()
        self.client_frontend[ip] = self.frontend.get_model()
        pass_filename = selrows[0]['invoce'].replace('/', '_')
        pass_data = self.print_pass(self.frontend.invoce_content_tab.options['rowData'])
        self.frontend.save_passport(pass_filename, pass_data)

        self.frontend.progressbar.visible = False

    async def on_send_oligos_button(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        selrows = await self.frontend.invoce_content_tab.get_selected_rows()

        if ip in list(self.oligomap_stack.input_selected_rows.keys()):
            self.oligomap_stack.input_selected_rows[ip].extend(selrows)
        else:
            self.oligomap_stack.input_selected_rows[ip] = []
            self.oligomap_stack.input_selected_rows[ip].extend(selrows)

        #print(len(self.oligomap_stack.input_selected_rows[ip]))

    def on_print_orders_date_range(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        start_date = datetime.strptime(self.frontend.start_date.value, "%Y-%m-%d")
        end_date = datetime.strptime(self.frontend.end_date.value, "%Y-%m-%d")

        data = pd.DataFrame(self.frontend.invoce_content_tab.options['rowData'])
        data['Date'] = pd.to_datetime(data['input date'], format='%m.%d.%Y')
        conditions = (data['Date'] >= start_date)&(data['Date'] <= end_date)
        columns = ["#", "Name", "5'-end", "Sequence", "3'-end", "Amount, oe", "Purification", "Lenght", "order id",
                   "client id", "input date"]
        df = data[conditions][columns]
        self.frontend.save_passport('data_range', df)

        self.client_frontend[ip] = self.frontend.get_model()


    def  print_pass(self, rowData):
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


    def load_history_ivoces(self):
        ip = self.app.storage.user.get('client_ip')
        self.pincode = self.client[ip]


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
            self.frontend.invoceChart.draw_stat_data(self.history_stat_data)

        self.client_frontend[ip] = self.frontend.get_model()


    def __getitem__(self, item):
        if item == 'on_load_button':
            return self.on_load_button
        elif item == 'on_show_actual_invoces_button':
            return self.on_show_actual_invoces_button
        elif item == 'on_pincode_change':
            return self.on_pincode_change
        elif item == 'on_show_invoces_content':
            return self.on_show_invoces_content
        elif item == 'on_show_by_status':
            return self.on_show_by_status
        elif item == 'on_print_invoce_passport':
            return self.on_print_invoce_passport
        elif item == 'on_send_oligos_button':
            return self.on_send_oligos_button
        elif item == 'on_print_orders_date_range':
            return self.on_print_orders_date_range
        elif item == 'on_update_invoces_tab':
            return self.on_update_invoces_tab