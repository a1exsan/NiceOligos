import pandas as pd
import requests
import json
from datetime import datetime

class api_db_interface():

    def __init__(self, db_IP, db_port):
        self.db_IP = db_IP
        self.db_port = db_port
        self.api_db_url = f'http://{self.db_IP}:{self.db_port}'
        self.pincode = ''
        self.client = {}
        self.client_frontend = {}

    def headers(self):
        return {'Authorization': f'Pincode {self.pincode}'}

class oligos_data_stack():
    def __init__(self):
        self.input_selected_rows = {}

class oligomaps_search(api_db_interface):

    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.pincode = ''
        self.maps_db_name = 'asm2000_map_1.db'
        self.strftime_format = "%Y-%m-%d"

    def map_in_progress(self, mapdata):
        map = json.loads(mapdata)
        df = pd.DataFrame(map)
        try:
            return round(df[df['Status'] == 'finished'].shape[0] * 100/df.shape[0], 0)
        except:
            return 100

    def get_oligomaps(self):
        self.all_not_fin_oligos = []
        url = f'{self.api_db_url}/get_all_tab_data/{self.maps_db_name}/main_map'
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            out = []
            for r in ret.json():
                d = {}
                d['#'] = r[0]
                d['Map name'] = r[2]
                d['Synth number'] = r[3]
                d['Date'] = r[1]
                d['in progress'] = self.map_in_progress(r[4])
                #d['map data'] = pd.DataFrame(json.loads(r[4]))

                df = pd.DataFrame(json.loads(r[4]))
                if 'Status' in list(df.keys()):
                    self.all_not_fin_oligos.extend(df[df['Status'] != 'finished'].to_dict('records'))
                    d['Finished'] = df[df['Status'] == 'finished'].shape[0]
                else:
                    d['Finished'] = df.shape[0]
                d['Total'] = df.shape[0]
                if 'Wasted' in list(df.keys()):
                    w_df = df[df['Wasted'] == True]
                    d['Wasted'] = w_df.shape[0]
                else:
                    d['Wasted'] = 0
                out.append(d)
            return out
        else:
            return []

    def get_oligomaps_data(self):
        url = f'{self.api_db_url}/get_all_tab_data/{self.maps_db_name}/main_map'
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            out = []
            for r in ret.json():
                d = {}
                d['#'] = r[0]
                d['Map name'] = r[2]
                d['Synth number'] = r[3]
                d['Date'] = r[1]
                d['in progress'] = self.map_in_progress(r[4])
                d['map data'] = pd.DataFrame(json.loads(r[4]))
                df = pd.DataFrame(d['map data'])
                if 'Status' in list(df.keys()):
                    d['Finished'] = df[df['Status'] == 'finished'].shape[0]
                else:
                    d['Finished'] = df.shape[0]
                d['Total'] = df.shape[0]
                if 'Wasted' in list(df.keys()):
                    d['Wasted'] = df[df['Wasted'] == True].shape[0]
                    #print(d['Wasted'])
                else:
                    d['Wasted'] = 0
                out.append(d)
            return out
        else:
            return []

    def get_actual_maps(self):
        total_maps = self.get_oligomaps_data()
        if len(total_maps) > 0:
            out = []
            for row in total_maps:
                df = row['map data']
                if df.shape[0] > 0:
                    if df[(df['DONE'] == True)|(df['Wasted'] == True)].shape[0] != df.shape[0]:
                        d = {}
                        d['#'] = row['#']
                        d['Map name'] = row['Map name']
                        d['Synth number'] = row['Synth number']
                        d['Date'] = row['Date']
                        d['in progress'] = row['in progress']
                        d['Finished'] = row['Finished']
                        d['Total'] = row['Total']
                        d['Wasted'] = row['Wasted']
                        out.append(d)
            return out
        else:
            return []

    def load_oligomap(self, seldata):
        if len(seldata) > 0:
            url = f"{self.api_db_url}/get_keys_data/{self.maps_db_name}/main_map/id/{seldata[0]['#']}"
            ret = requests.get(url, headers=self.headers())
            if ret.status_code == 200:
                self.oligo_map_id = seldata[0]['#']
                meta = ret.json()
                map_data = json.loads(meta[0][4])
                accord_data = json.loads(meta[0][5])
                #print(meta)
                return map_data, accord_data, meta[0][2], meta[0][3]
            else:
                return [], []

    def get_oligomaps_data_tail(self, tail_len):
        url = f'{self.api_db_url}/get_all_tab_data_tail/{self.maps_db_name}/main_map/{tail_len}'
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            out = []
            for r in ret.json():
                d = {}
                d['#'] = r[0]
                d['Map name'] = r[2]
                d['Synth number'] = r[3]
                d['Date'] = r[1]
                d['map data'] = json.loads(r[4])
                out.append(d)
            return out
        else:
            return []

    def get_oligomaps_date_tail(self, date, tail_len=20):
        url = f'{self.api_db_url}/get_oligomaps_date_tail/{self.maps_db_name}/main_map/{date}/{tail_len}'
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            out = []
            for r in ret.json():
                d = {}
                d['#'] = r['0']
                d['Map name'] = r['2']
                d['Synth number'] = r['3']
                d['Date'] = r['1']
                d['map data'] = json.loads(r['4'])
                out.append(d)
            return out
        else:
            return []

    def find_amount_by_order_id(self, order_id):
        maps = pd.DataFrame()
        for map in self.map_list['map data']:
            maps = pd.concat([maps, pd.DataFrame(map)])
        maps.reset_index(inplace=True)
        maps = maps[(maps['Order id'] == int(order_id))&(maps['Wasted'] == False)]
        return maps

    def get_low_amount_limit(self, amount):
        if amount.find('-')>-1:
            return float(amount[:amount.find('-')])
        else:
            return float(amount)

    def insert_map_to_base(self, name, synth_number, rowData, accordData):
        if len(rowData) > 0:

            synth_date = datetime.now().date().strftime(self.strftime_format)

            url = f"{self.api_db_url}/insert_data/{self.maps_db_name}/main_map"
            r = requests.post(url,
                              json=json.dumps([synth_date, name, synth_number,
                                               json.dumps(rowData),
                                               json.dumps(accordData)]), headers=self.headers())

            url = f'{self.api_db_url}/get_all_tab_data/{self.maps_db_name}/main_map'
            ret = requests.get(url, headers=self.headers())
            indexs = [r[0] for r in ret.json()]

            self.oligo_map_id = max(indexs)

            data = []
            for row in rowData:
                d = row.copy()
                d['map #'] = self.oligo_map_id
                data.append(d)

            url = f"{self.api_db_url}/update_data/{self.maps_db_name}/main_map/{self.oligo_map_id}"
            r = requests.put(url,
                             json=json.dumps({
                                 'name_list': ['map_tab'],
                                 'value_list': [
                                     json.dumps(data)
                                 ]
                             })
                             , headers=self.headers())

            return r.status_code
        else:
            return 404


    def get_order_status(self, row):
        state_list = ['synth', 'sed', 'click', 'cart', 'hplc', 'paag', 'LCMS', 'subl']
        flag_list = []
        for state in state_list:
            flag_list.append(row[f'Do {state}'] == row[f'Done {state}'])
        status = 'synthesis'
        for i in range(8):
            if not flag_list[i]:
                if i < 3:
                    status = 'synthesis'
                elif i > 2 and i < 6:
                    status = 'purification'
                elif i == 7:
                    status = 'formulation'
                return status
            else:
                if i == 7:
                    status = 'finished'
                    return status
        return status


    def update_oligomap_status(self, rowData, accordrowdata):
        if len(rowData) > 0:
            print(rowData[0])
            if 'map #' in list(rowData[0].keys()):
                for row in rowData:
                    if row['map #'] != '':
                        self.oligo_map_id = int(row['map #'])
                        print(f'MAP ID: {self.oligo_map_id}')
                        break
        if self.oligo_map_id > -1:
            out = []
            for row in rowData:
                out.append(row)
                out[-1]['Date'] = datetime.now().date().strftime('%d.%m.%Y')
                out[-1]['Status'] = self.get_order_status(row)
                if out[-1]['Status'] == 'finished':
                    out[-1]['DONE'] = True
                else:
                    out[-1]['DONE'] = False

            url = f"{self.api_db_url}/update_data/{self.maps_db_name}/main_map/{self.oligo_map_id}"
            r = requests.put(url,
                              json=json.dumps({
                                  'name_list': ['map_tab', 'accord_tab'],
                                  'value_list': [
                                      json.dumps(out),
                                      json.dumps(accordrowdata)
                                  ]
                              })
                             , headers=self.headers())
            print(f'update status {self.oligo_map_id}: {r.status_code}')
            return out
        else:
            return rowData



def test():
    osearch = oligomaps_search('127.0.0.1', '8012')
    osearch.pincode = '4720'
    #map_list = osearch.get_oligomaps_data_tail(40)
    #for row in map_list:
    #    print(row)

    map_list = osearch.get_oligomaps_date_tail('2025-07-02', 20)
    print(pd.DataFrame(map_list))

def test2():
    osearch = oligomaps_search('127.0.0.1', '8012')
    osearch.pincode = '4720'
    maps = osearch.find_amount_by_order_id(6205, '2025-07-02')
    print(maps.keys())
    df = maps['Dens, oe/ml'] * maps['Vol, ml']
    print(df.sum())

if __name__ == '__main__':
    test2()