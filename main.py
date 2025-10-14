from nicegui import ui, app
from typing import Optional
from fastapi.responses import RedirectResponse
from OligoMap_utils import oligos_data_stack

from invoce_page import invoice_page_model
from invoce_page import navigation_menu
from oligosynth_page import oligosynth_panel_page_model
from raw_material_page import rawmaterial_panel_page_model
from chemicals_page import chemicals_page_model
from input_order_page import input_order_page_model
from invoce_page import api_db_interface
from stock_data_page import stock_data_page_model
from molseq_lang import modification_page_model
from server import IP

import hashlib

oligo_map_stack = oligos_data_stack()

users = {'test': hashlib.sha256('12345'.encode('utf-8')).hexdigest()}

#IP_addr = '127.0.0.1'
#IP_addr = '192.168.16.145'
#IP_addr = '84.252.133.233'
IP_addr = IP

app.storage.general['db_IP'] = IP_addr
app.storage.general['db_port'] = '8012'
api_object = api_db_interface(app.storage.general.get('db_IP'), app.storage.general['db_port'])
app.add_static_files('/img', 'static_images')


@ui.page('/')
def main_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    with ui.column():
        navi_front = navigation_menu(IP_addr, '8012')
        ui.button('Выйти', on_click=logout, icon='logout').props('outline round')
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: 1000px')
        #ui.label(f'Добро пожаловать, {app.storage.user["email"]}!').classes('text-2xl')

def logout():
    app.storage.user.clear()
    ui.link('Войти', '/login')
    ui.run_javascript(
        '''
        document.querySelector('a[href="/login"]').click();
        '''
    )

@ui.page('/register')
def register_page():
    ui.dark_mode(True)
    email = ui.input('Email')
    password = ui.input('Пароль', password=True)
    confirm_password = ui.input('Подтверждение пароля', password=True)

    def register():
        if not email.value or not password.value:
            ui.notify('Поля не должны быть пустыми', color='negative')
            return
        if password.value != confirm_password.value:
            ui.notify('Пароли не совпадают', color='negative')
            return
        if email.value in users:
            ui.notify('Пользователь с таким email уже зарегистрирован', color='negative')
            return
        users[email.value] = hashlib.sha256(password.value.encode('utf-8')).hexdigest()
        ui.notify('Регистрация успешна!', color='positive')
        ui.link('Войти', '/login')

    ui.button('Зарегистрироваться', on_click=register)
    ui.link('Есть аккаунт? Войти', '/login')

@ui.page('/login')
def login_page() -> Optional[RedirectResponse]:
    if app.storage.user.get('authenticated', False):
        return RedirectResponse('/')
    ui.dark_mode(True)

    def login():

        ps = f'{email.value}:{password.value}'
        api_object.pincode = hashlib.sha256(ps.encode('utf-8')).hexdigest()
        #print(api_object.check_pincode())
        #if users.get(email.value) == hashlib.sha256(password.value.encode('utf-8')).hexdigest():
        if api_object.check_pincode().status_code == 200:
            app.storage.user.update({'email': email.value, 'authenticated': True})
            app.storage.user.update({'pincode': api_object.pincode})
            ret = api_object.get_user_data()
            user_data = ret.json()
            app.storage.user.update({'user_status': user_data['status']})
            ui.link('Главная', '/')
            ui.run_javascript(
                '''
                document.querySelector('a[href="/"]').click();
                '''
            )
        else:
            ui.notify('Неверный email или пароль', color='negative')

    with ui.row():
        email = ui.input('Email')
        password = ui.input('Пароль', password=True)
    ui.button('Войти', on_click=login)
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: 1000px;')


#@ui.page('/main')
#def index():
#    navi_front = navigation_menu(IP_addr, '8012')
#    ui.image('images/background_1.jpeg').style('max-width: 100%; height: auto;')

@ui.page('/input_order_panel')
def input_order_panel() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner', 'synth_master']:
        order_page = input_order_page_model(IP_addr, '8012')

@ui.page('/invoce_panel')
def invoce_panel_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner', 'lab_master', 'synth_master']:
        invoce_front = invoice_page_model(IP_addr, '8012')

@ui.page('/oligosynth_panel')
def oligosynth_panel_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner', 'lab_master', 'synth_master']:
        oligosynt_front = oligosynth_panel_page_model(IP_addr, '8012')


@ui.page('/rawmaterials_panel')
def rawmaterials_panel_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner', 'lab_master', 'synth_master']:
        rawmat_panel = rawmaterial_panel_page_model(IP_addr, 8012)

@ui.page('/chemicals_panel')
def chemicals_panel_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner']:
        chem_page = chemicals_page_model(IP_addr, 8012)

@ui.page('/stock_data')
def stock_data_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner', 'stock_viewer']:
        model = stock_data_page_model()


@ui.page('/modifications')
def main_page() -> Optional[RedirectResponse]:
    if not app.storage.user.get('authenticated', False):
        return RedirectResponse('/login')#

    navi_front = navigation_menu(IP_addr, '8012')
    if app.storage.user.get('user_status') in ['own', 'owner']:
       model = modification_page_model()


#if __name__ == '__main__':
ui.run(
        storage_secret='lfjvnieurvqneorivnq/lke',
        title='NiceOligo',
        favicon="static_images/favicon.ico",
        port=8080
    )
