import json

from nicegui import ui, app
from OligoMap_utils import api_db_interface
import requests
import hashlib
from datetime import datetime

class dialog_confirm():
    def __init__(self):
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                with ui.row():
                    ui.label('Вы уверены !!!')
                    self.ch_btn = ui.button('Да', on_click=self.on_enter)
                    ui.button('Нет', on_click=self.dialog.close)

    def return_func(self):
        pass

    def on_enter(self):
        self.return_func()
        self.dialog.close()

class dialog_add_user():
    def __init__(self):
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                with ui.row():
                    self.login = ui.input('Login:', password=False).style('font-size: 24px;')
                    self.pwd = ui.input('pass:', password=False).style('font-size: 24px;')
                    self.ch_btn = ui.button('ввод', on_click=self.on_enter)
                    self.ch_btn.on('keydown.enter', self.on_enter)


    def return_func(self, login, password):
        pass

    def on_enter(self):
        self.return_func(self.login.value, self.pwd.value)
        self.dialog.close()

class user_admin_model(api_db_interface):
    def __init__(self):
        self.db_IP = app.storage.general.get('db_IP')
        self.db_port = app.storage.general.get('db_port')
        super().__init__(self.db_IP, self.db_port)
        self.pincode = app.storage.user.get('pincode')
        self.users_tab = None
        self.new_user_pincode = None
        self.new_user_login = None
        self.sel_user_id = None
        self.pass_dialog = dialog_add_user()
        self.pass_dialog.return_func = self.pass_func
        self.confirm_dialog = dialog_confirm()
        self.confirm_dialog.return_func = self.delete_user_row

        self.users_rowdata = self.get_users_data()
        with ui.column():
            with ui.row():
                ui.button('show users', on_click=self.show_users_tab)
                ui.button('update', on_click=self.update_users_data)
                ui.button('add user', color='green', on_click=self.add_user)
                ui.button('delete user', color='red', on_click=self.delete_user)
            self.get_tab(self.users_rowdata)


    def get_tab(self, rowdata):
        columnDefs = [
            {"field": "id", 'editable': False},
            {"field": "login", 'editable': True},
            {"field": "pass", 'editable': False},
            {"field": "pass_date", 'editable': True},
            {"field": "data_json", 'editable': True},
        ]

        self.users_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': rowdata,
                'rowSelection': 'multiple',
                "pagination": True,
                "enterNavigatesVertically": True,
                "enterNavigatesVerticallyAfterEdit": True,
                "singleClickEdit": True
                # ':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                #                ' { background: "green" }',
            }
            ,
            theme='alpine-dark').style("width: 1200px; height: 800px")  # alpine  material  quartz  balham
        self.users_tab.auto_size_columns = True
        self.users_tab.on("cellValueChanged", self.users_tab_update_cell)

    def get_users_data(self):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_users}/users'
        ret = requests.get(url, headers=self.headers())
        rowdata = []
        for row in ret.json():
            d = {}
            d['id'] = row[0]
            d['login'] = row[1]
            d['pass'] = row[2]
            d['pass_date'] = row[3]
            d['data_json'] = row[4]
            rowdata.append(d)
        return rowdata

    def update_users_data_json(self, id, data_json):
        url = f'{self.api_db_url}/update_data/{self.db_users}/users/{id}'
        ret = requests.put(url,
                           json=json.dumps({
                               'name_list': ['data_json'],
                               'value_list': [data_json]
                           }),
                           headers=self.headers())
        if ret.status_code == 200:
            ui.notify('Данные успешно обновлены')
        else:
            ui.notify('Данные не удалось обновить')

    def users_tab_update_cell(self, e):
        self.users_tab.options['rowData'][e.args["rowIndex"]] = e.args["data"]

    def show_users_tab(self):
        self.users_rowdata = self.get_users_data()
        self.users_tab.options['rowData'] = self.users_rowdata
        self.users_tab.update()

    def update_users_data(self):
        for row in self.users_tab.options['rowData']:
            self.update_users_data_json(row['id'], row['data_json'])
        self.users_rowdata = self.get_users_data()
        self.users_tab.options['rowData'] = self.users_rowdata
        self.users_tab.update()

    def pass_func(self, login, pwd):
        s = f'{login}:{pwd}'
        self.new_user_pincode = hashlib.sha256(s.encode('utf-8')).hexdigest()
        self.new_user_login = login

        url = f'{self.api_db_url}/add_new_user/{self.db_users}/users'
        ret = requests.post(url,
                            json=json.dumps([
                                login, self.new_user_pincode, datetime.now().date().strftime('%d.%m.%Y'),
                                json.dumps({'status': 'synth_master'})
                            ]),
                            headers=self.headers())
        if ret.status_code == 200:
            ui.notify('Пользователь добавлен')
            self.show_users_tab()

    def add_user(self):
        self.pass_dialog.dialog.open()

    async def delete_user(self):
        selrows = await self.users_tab.get_selected_rows()
        if len(selrows) > 0:
            self.sel_user_id = selrows[0]['id']
            self.confirm_dialog.dialog.open()

    def delete_user_row(self):
        if self.sel_user_id is not None:
            url = f'{self.api_db_url}/delete_data/{self.db_users}/users/{self.sel_user_id}'
            ret = requests.delete(url, headers=self.headers())
            if ret.status_code == 200:
                ui.notify('Данные успешно удалены')
                self.show_users_tab()