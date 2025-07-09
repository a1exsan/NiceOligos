import pandas as pd


class menu_front():
    def __init__(self, ui):
        with ui.button(icon='menu'):
            with ui.menu():
                ui.menu_item('Invoces', on_click=lambda: ui.open('https://example.com'))
                ui.menu_item('Option 2')
                with ui.menu_item('Option 3', auto_close=False):
                    with ui.item_section().props('side'):
                        ui.icon('keyboard_arrow_right')
                    with ui.menu().props('anchor="top end" self="top start" auto-close'):
                        ui.menu_item('Sub-option 1')
                        ui.menu_item('Sub-option 2')
                        ui.menu_item('Sub-option 3')

class navigation_menu():
    def __init__(self, ui):
        dark = ui.dark_mode(True)
        #ui.switch('Dark mode').bind_value(dark)

        navigation = ui.row()
        with navigation:
            ui.link('Home', '/').style('font-size: 24px;')
            ui.link('Invoces', '/invoce_panel').style('font-size: 24px;')


class invoice_frontend():
    def __init__(self, ui):

        rowData = pd.DataFrame({
            '#': [0],
            'invoce': [''],
            'client': [''],
            'input date': [''],
            'out date': [''],
            'number': [0],
            'in queue%': [''],
            'synth%': [''],
            'purif%': [''],
            'formul%': [''],
            'fin%': [''],
            'product days': ['3'],
            'status': [''],
            'send': [True],
            'value P': ['']
        })

        colDefs = [
                        {"field": "#", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "invoce", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "client", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "input date", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "out date", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "number", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "in queue%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "synth%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "purif%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "formul%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "fin%", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "product days", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "status", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "send", 'editable': True, 'filter': 'agSetColumnFilter'},
                        {"field": "value P", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                    ]

        self.input_pincode = ui.input(label='pincode', placeholder='enter pincode')

        summary_panel = ui.row()
        with summary_panel:
            col1 = ui.column()
            with col1:
                self.total_oligos_invoces = ui.circular_progress(min=0, max=1000, value=445,
                                                                 size='96px', color='orange')
                ui.label('Total oligos')
            col2 = ui.column()
            with col2:
                self.in_queue_oligos_invoces = ui.circular_progress(min=0, max=1000, value=100,
                                                                    size='96px', color='orange')
                ui.label('In queue')
            col3 = ui.column()
            with col3:
                self.synth_oligos_invoces = ui.circular_progress(min=0, max=1000, value=45,
                                                                 size='96px', color='orange')
                ui.label('Synthesis')
            col4 = ui.column()
            with col4:
                self.purif_oligos_invoces = ui.circular_progress(min=0, max=1000, value=65,
                                                                 size='96px', color='orange')
                ui.label('Purification')
            col5 = ui.column()
            with col5:
                self.formul_oligos_invoces = ui.circular_progress(min=0, max=1000, value=180,
                                                                  size='96px', color='orange')
                ui.label('Formulation')
            col6 = ui.column()
            with col6:
                self.fin_oligos_invoces = ui.circular_progress(min=0, max=1000, value=230,
                                                               size='96px', color='orange')
                ui.label('Finished')
            col7 = ui.column()
            with col7:
                self.total_price_oligos_invoces = ui.circular_progress(min=0, max=5000000, value=880000,
                                                                       size='96px', color='orange')
                ui.label('Price')

        self.ag_grid = ui.aggrid(
            {
                'columnDefs': colDefs,
                'rowData': rowData.to_dict('records'),
                'rowSelection': 'multiple'
            }
        ,
        theme='alpine-dark').classes('h-[800px]') # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True

        with ui.row():
            self.on_load_button = ui.button('load invoces')
            self.on_show_actual_invoces_button = ui.button('Show actual')

    def __getitem__(self, item):
        if item == 'on_load_button':
            return self.on_load_button
        elif item == 'on_show_actual_invoces_button':
            return self.on_show_actual_invoces_button

    def get_element_list(self, key):
        if key == 'button':
            return [
                'on_load_button',
                'on_show_actual_invoces_button'
            ]
        else:
            return []