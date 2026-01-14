from nicegui import app, ui
from OligoMap_utils import api_db_interface
import hashlib

class dialog_text_input():
    def __init__(self):
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                with ui.row():
                    self.pwd = ui.input('Подтвердите старым:', password=True).style('font-size: 24px;')
                    self.ch_btn = ui.button('ввод', on_click=self.on_enter)
                    self.ch_btn.on('keydown.enter', self.on_enter)


    def return_func(self, text):
        pass

    def on_enter(self):
        self.return_func(self.pwd.value)
        self.dialog.close()


class user_admin_profile():
    def __init__(self):
        self.check_dialog = dialog_text_input()
        self.check_dialog.return_func = self.check_pwd

        with ui.column():
            with ui.row():
                self.new_pwd = ui.input(label='Новый пароль:', password=True)
                self.new_pwd.props('autocomplete="new-password"')
                self.repeat_pwd = ui.input(label='Повторить новый пароль:', password=True)
                self.repeat_pwd.props('autocomplete="new-password"')
                self.repeat_pwd.on('keydown.enter', self.change_pwd)
            ui.button(text='Сохранить новый пароль', on_click=self.change_pwd)

    def check_pwd(self, text):
        user_email = app.storage.user.get("email")
        ps = f'{user_email}:{text}'
        pincode = hashlib.sha256(ps.encode('utf-8')).hexdigest()
        if pincode == app.storage.user.get('pincode'):
            dbi = api_db_interface(app.storage.general.get('db_IP'), app.storage.general['db_port'])
            dbi.pincode = pincode

            ps_new = f'{user_email}:{self.new_pwd.value}'
            pincode_new = hashlib.sha256(ps_new.encode('utf-8')).hexdigest()

            ret = dbi.change_pwd(user_email, pincode, pincode_new)
            #ui.notify(f'{ret}')
            if ret.status_code == 200:
                ui.notify('Пароль успешно изменен')
                app.storage.user.clear()
                ui.run_javascript('window.location.href = "/login"')


    def change_pwd(self):
        if self.new_pwd.value == self.repeat_pwd.value and self.new_pwd.value != '':
            self.check_dialog.dialog.open()
        elif self.new_pwd.value == '':
            ui.notify('Пароль не должен быть пустым')
        else:
            ui.notify('Пароли не совпадают')



