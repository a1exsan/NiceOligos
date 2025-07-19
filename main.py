from nicegui import ui, app
import pandas as pd

import frontend
import invoce_backend

invoce_tab = invoce_backend.invoce_table('127.0.0.1', '8012')
#invoce_tab = invoce_backend.invoce_table('192.168.16.145', '8012')


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

@ui.page('/oligosynth_panel')
async def oligosynth_panel_page(client):
    await client.connected()
    ip = client.environ['asgi.scope']['client']
    app.storage.user['client_ip'] = ip[0]
    navi_front = frontend.navigation_menu(ui)


ui.run(storage_secret='NiceGUI_oligo_app_1')
