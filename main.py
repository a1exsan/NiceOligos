from nicegui import ui, app
from OligoMap_utils import oligos_data_stack

from frontend import invoice_page_model
from frontend import navigation_menu
from frontend_oligosynth_panel import oligosynth_panel_page_model
from raw_mat_frontend import rawmaterial_panel_page_model

oligo_map_stack = oligos_data_stack()

IP_addr = '127.0.0.1'
#IP_addr = '192.168.16.145'

app.add_static_files('/img', 'static_images')

@ui.page('/')
def index():
    navi_front = navigation_menu(IP_addr, '8012')
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: auto;')

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


ui.run(storage_secret='NiceGUI_oligo_app_1')
