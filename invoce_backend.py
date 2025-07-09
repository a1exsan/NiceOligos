import pandas as pd
import requests

from datetime import datetime
from datetime import timedelta

class api_db_interface():

    def __init__(self, db_IP, db_port):
        self.db_IP = db_IP
        self.db_port = db_port
        self.api_db_url = f'http://{self.db_IP}:{self.db_port}'
        self.pincode = ''

    def headers(self):
        return {'Authorization': f'Pincode {self.pincode}'}

class invoce_table(api_db_interface):

    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.pincode = ''
        self.db_name = 'scheduler_oligolab_2.db'


    def init_frontend(self, front):
        self.frontend = front
        self.frontend.input_pincode.value = self.pincode


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


    def on_load_button(self):
        self.pincode = self.frontend.input_pincode.value

        self.frontend.ag_grid.options['rowData'] = self.get_all_invoces()
        self.frontend.ag_grid.update()


    def on_show_actual_invoces_button(self):
        self.pincode = self.frontend.input_pincode.value

        df = pd.DataFrame(self.get_all_invoces())
        df = df[df['status'] == 'in progress']

        self.frontend.ag_grid.options['rowData'] = self.set_invoces_timing_prognosis(df.to_dict('records'))
        self.frontend.ag_grid.update()

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


    def __getitem__(self, item):
        if item == 'on_load_button':
            return self.on_load_button
        elif item == 'on_show_actual_invoces_button':
            return self.on_show_actual_invoces_button