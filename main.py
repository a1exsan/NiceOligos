from nicegui import ui, app
from OligoMap_utils import oligos_data_stack

from invoce_page import invoice_page_model
from invoce_page import navigation_menu
from oligosynth_page import oligosynth_panel_page_model
from raw_material_page import rawmaterial_panel_page_model
from chemicals_page import chemicals_page_model
from input_order_page import input_order_page_model

oligo_map_stack = oligos_data_stack()

#IP_addr = '127.0.0.1'
IP_addr = '192.168.16.145'

app.add_static_files('/img', 'static_images')

@ui.page('/')
def index():
    navi_front = navigation_menu(IP_addr, '8012')
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: auto;')

@ui.page('/input_order_panel')
def input_order_panel():
    navi_front = navigation_menu(IP_addr, '8012')
    order_page = input_order_page_model(IP_addr, '8012')

@ui.page('/invoce_panel')
def invoce_panel_page():
    navi_front = navigation_menu(IP_addr, '8012')
    invoce_front = invoice_page_model(IP_addr, '8012')

@ui.page('/oligosynth_panel')
def oligosynth_panel_page(client):
    navi_front = navigation_menu(IP_addr, '8012')
    oligosynt_front = oligosynth_panel_page_model(IP_addr, '8012')


@ui.page('/rawmaterials_panel')
def rawmaterials_panel_page():
        navi_front = navigation_menu(IP_addr, '8012')
        rawmat_panel = rawmaterial_panel_page_model(IP_addr, 8012)

@ui.page('/chemicals_panel')
def chemicals_panel_page():
        navi_front = navigation_menu(IP_addr, '8012')
        chem_page = chemicals_page_model(IP_addr, 8012)


ui.run(storage_secret='NiceGUI_oligo_app_1')
