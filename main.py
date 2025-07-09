from nicegui import ui
import pandas as pd

import frontend
import invoce_backend

invoce_tab = invoce_backend.invoce_table('127.0.0.1', '8012')
#invoce_tab = invoce_backend.invoce_table('192.168.16.145', '8012')


@ui.page('/')
def index():
    navi_front = frontend.navigation_menu(ui)
    ui.image('images/background_1.jpeg').style('max-width: 100%; height: auto;')

@ui.page('/invoce_panel')
def invoce_panel_page():
    navi_front = frontend.navigation_menu(ui)
    invoce_front = frontend.invoice_frontend(ui)
    invoce_tab.init_frontend(invoce_front)
    # ADD Invoce page buttons
    for btn in invoce_front.get_element_list('button'):
        invoce_front[btn].on('click', invoce_tab[btn])

ui.run()
