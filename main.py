from nicegui import ui, app
from OligoMap_utils import oligos_data_stack

import frontend
import invoce_backend
import frontend_oligosynth_panel
import Synthesis_builder_backend
import raw_mat_frontend
import stock_backend

oligo_map_stack = oligos_data_stack()

IP_addr = '127.0.0.1'
#IP_addr = '192.168.16.145'

invoce_tab = invoce_backend.invoce_table(IP_addr, '8012', oligo_map_stack)
oligomap_tab = Synthesis_builder_backend.Oligomap_backend(IP_addr, '8012', oligo_map_stack)
stock_tab = stock_backend.stock_backend_model(IP_addr, '8012')


app.add_static_files('/img', 'static_images')

@ui.page('/')
async def index(client):
    await client.connected()
    ip = client.environ['asgi.scope']['client']
    app.storage.user['client_ip'] = ip[0]
    invoce_tab.init_application(app)
    print(ip, ip[0])
    navi_front = frontend.navigation_menu(ui)
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: auto;')

@ui.page('/invoce_panel')
async def invoce_panel_page(client):
    await client.connected()
    ip = client.environ['asgi.scope']['client']
    #print(ip)
    navi_front = frontend.navigation_menu(ui)
    invoce_tab.init_application(app)
    invoce_front = frontend.invoice_frontend(ui)
    invoce_tab.init_frontend(invoce_front, ip[0])
    # ADD Invoce page buttons
    for btn, func in invoce_front.get_element_list('button'):
        invoce_front[btn].on(func, invoce_tab[btn])
    invoce_front['on_pincode_change'].on_value_change(invoce_tab['on_pincode_change'])
    invoce_front['on_show_by_status'].on_value_change(invoce_tab['on_show_by_status'])

    #print('INVOCE', oligomap_tab.client_frontend[ip[]])

@ui.page('/oligosynth_panel')
async def oligosynth_panel_page(client):
    await client.connected()
    ip = client.environ['asgi.scope']['client']
    app.storage.user['client_ip'] = ip[0]
    navi_front = frontend.navigation_menu(ui)

    oligosynt_front = frontend_oligosynth_panel.oligosynth_panel(oligomap_tab)
    oligomap_tab.init_frontend(oligosynt_front, ip[0], invoce_tab.get_pincode(ip[0]))

    for btn, func in oligosynt_front.get_element_list('button'):
        oligosynt_front[btn].on(func, oligomap_tab[btn])
    oligosynt_front['synth_scale_selector'].on_value_change(oligomap_tab['synth_scale_selector'])
    oligosynt_front['wells_layer_selector'].on_value_change(oligomap_tab['wells_layer_selector'])
    oligosynt_front.done_event = oligomap_tab.on_sel_done_btn
    oligosynt_front.do_event = oligomap_tab.on_sel_do_btn

@ui.page('/rawmaterials_panel')
async def rawmaterials_panel_page(client):
        await client.connected()
        ip = client.environ['asgi.scope']['client']
        navi_front = frontend.navigation_menu(ui)

        rawmat_panel = raw_mat_frontend.rawmaterial_panel(stock_tab, IP_addr, 8012)
        rawmat_panel.pincode = invoce_tab.get_pincode(ip[0])
        rawmat_panel.init_frontend()


ui.run(storage_secret='NiceGUI_oligo_app_1')
