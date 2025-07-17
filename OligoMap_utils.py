import pandas as pd
import requests
import json

class api_db_interface():

    def __init__(self, db_IP, db_port):
        self.db_IP = db_IP
        self.db_port = db_port
        self.api_db_url = f'http://{self.db_IP}:{self.db_port}'
        self.pincode = ''

    def headers(self):
        return {'Authorization': f'Pincode {self.pincode}'}


class oligomaps_search(api_db_interface):

    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.pincode = ''
        self.maps_db_name = 'asm2000_map_1.db'

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