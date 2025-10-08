from OligoMap_utils import api_db_interface
import requests
from nicegui import ui, app
import json

class synth_base(api_db_interface):
    def __init__(self):
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        super().__init__(IP, port)
        self.pincode = app.storage.user.get('pincode')
        self.selected_id = 0

    def get_all_data(self):
        url = f'{self.api_db_url}/get_all_tab_data/{self.synth_db_name}/main'
        r = requests.get(url, headers=self.headers())
        if r.status_code == 200:
            return r.json()
        else:
            return []

    def get_all_rowdata(self):
        out = []
        tab = self.get_all_data()
        for row in tab:
            d = {
                'id': row[0],
                'synth_name': row[1],
                'scale': row[2],
                'data_json': row[3],
                'template': row[4]
            }
            out.append(d)
        return out

    def insert_method(self, data_object):
        tab = self.get_all_data()
        insert = True
        for row in tab:
            if row[3] == data_object['data_json']:
                insert = False
                break
        if insert:
            url = f'{self.api_db_url}/insert_data/{self.synth_db_name}/main'
            r = requests.post(url,
                          json=json.dumps([
                              data_object['synth_name'],
                              data_object['scale'],
                              data_object['data_json'],
                              data_object['template']
                          ]), headers=self.headers())
            if r.status_code == 200:
                ui.notify(f'Метод добавлен в базу')
        else:
            ui.notify(f'Метод есть в базе')

    def update_method(self, id, data):
        print(id, data)
        url = f'{self.api_db_url}/update_data/{self.synth_db_name}/main/{id}'
        r = requests.put(url,
                         json=json.dumps({
                             'name_list': list(data.keys()),
                             'value_list': list(data.values())
                         }), headers=self.headers())
        if r.status_code == 200:
            ui.notify(f'Метод успешно обновлен')