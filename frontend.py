import pandas as pd
import matplotlib.pyplot as plt

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
            ui.link('Oligo synthesis', '/oligosynth_panel').style('font-size: 24px;')


class invoice_frontend():
    def __init__(self, ui):

        self.ui = ui
        self.model = {}

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

        self.on_pincode_change = ui.input(label='pincode', placeholder='enter pincode')

        summary_panel = ui.row()
        with (summary_panel):
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
            #self.invoces_stat_mage = ui.image()
            #self.invoces_stat_mage.set_source('images/number_oligos_plot_1.png')
            #self.invoces_stat_mage.props('width=1250px height=180px')

            #plt.style.use('dark_background')
            #with ui.matplotlib(figsize=(15, 4)).figure as self.invoce_stat_figure:
            #    x = [1,2,3]
            #    y = [1,2,3]
            #    self.invoce_stat_ax = self.invoce_stat_figure.gca()
            #    self.invoce_stat_ax.hist(x, y)
            #    self.invoce_stat_ax.set_xlabel('Month')
            #    self.invoce_stat_ax.set_ylabel('Number of oligos')


        self.ag_grid = ui.aggrid(
            {
                'columnDefs': colDefs,
                'rowData': rowData.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                #"enableRangeSelection": True,
            }
        ,
        theme='alpine-dark').classes('h-[800px]') # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True

        with ui.row():
            self.on_load_button = ui.button('load invoces')
            self.on_show_actual_invoces_button = ui.button('Show actual')
            self.on_show_invoces_content = ui.button('invoce content')
            self.on_print_invoce_passport = ui.button('print passport', color='#FFA000')
            self.progressbar = ui.spinner(size='md', color='#FFA000')
        self.progressbar.visible = False
        self.on_load_button.props('id="on_load_button"')

        status_list = [
            'in queue',
            'synthesis',
            'purification',
            'formulation',
            'finished',
            'in progress',
            'wasted in progress',
            'not compleated in finished',
        ]

        with ui.row():
            self.on_show_by_status = ui.select(options=status_list, with_input=True).classes('w-60')

        invoce_content_df = pd.DataFrame(
            {
                '#': [1],
                'Name': [''],
                "5'-end": [''],
                'Sequence': [''],
                "3'-end": [''],
                'Amount, oe': ['5-10'],
                'Exist, oe': [0.],
                'Purification': ['Cart'],
                'Lenght': [''],
                'status': ['in queue'],
                'input date': [''],
                'output date': [''],
                'client id': [''],
                'order id': [''],
                'sufficiency': [0.]
            }
        )

        columnDefs = [
            {
                "field": "#",
                "checkboxSelection": True,
                "headerCheckboxSelection": True,
                "headerCheckboxSelectionFilteredOnly": True,
            },
            {"field": "Name", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "5'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Sequence", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "3'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Amount, oe", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Exist, oe", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Purification", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Lenght", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "status", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "input date", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "output date", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "order id", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "client id", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "sufficiency", 'filter': 'agTextColumnFilter', 'floatingFilter': True}
        ]

        self.invoce_content_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': invoce_content_df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                #"enterNavigatesVertically": True,
                #"enterNavigatesVerticallyAfterEdit": True,
                #"singleClickEdit": True
                ':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                                ' { background: "green" }',
            }
        ,
        theme='alpine-dark').classes('h-[800px]') # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True

    def __getitem__(self, item):
        if item == 'on_load_button':
            return self.on_load_button
        elif item == 'on_show_actual_invoces_button':
            return self.on_show_actual_invoces_button
        elif item == 'on_pincode_change':
            return self.on_pincode_change
        elif item == 'on_show_invoces_content':
            return self.on_show_invoces_content
        elif item == 'on_show_by_status':
            return self.on_show_by_status
        elif item == 'on_print_invoce_passport':
            return self.on_print_invoce_passport

    def get_element_list(self, key):
        if key == 'button':
            return [
                ('on_load_button', 'click'),
                ('on_show_actual_invoces_button', 'click'),
                ('on_show_invoces_content', 'click'),
                ('on_print_invoce_passport', 'click'),
            ]
        else:
            return []

    def save_passport(self, filename, data):
        csv_data = data.to_csv(index=False, sep='\t').encode('utf-8')
        #csv_data = data.to_excel(sheet_name='Sheet1', index=False, sep='\t').encode('utf-8')
        self.ui.download(csv_data, f'{filename}.csv')
        #self.ui.download(csv_data, f'{filename}.xslx')


    def get_model(self):
        self.model['on_pincode_change'] = self.on_pincode_change
        self.model['total_oligos_invoces'] = self.total_oligos_invoces
        self.model['in_queue_oligos_invoces'] = self.in_queue_oligos_invoces
        self.model['synth_oligos_invoces'] = self.synth_oligos_invoces
        self.model['purif_oligos_invoces'] = self.purif_oligos_invoces
        self.model['formul_oligos_invoces'] = self.formul_oligos_invoces
        self.model['fin_oligos_invoces'] = self.fin_oligos_invoces
        self.model['total_price_oligos_invoces'] = self.total_price_oligos_invoces
        self.model['ag_grid'] = self.ag_grid
        self.model['invoce_content_tab'] = self.invoce_content_tab
        return self.model

    def set_model(self, model):
        self.on_pincode_change.value = model['on_pincode_change'].value
        self.total_oligos_invoces.value = model['total_oligos_invoces'].value
        self.in_queue_oligos_invoces.value = model['in_queue_oligos_invoces'].value
        self.synth_oligos_invoces.value = model['synth_oligos_invoces'].value
        self.purif_oligos_invoces.value = model['purif_oligos_invoces'].value
        self.formul_oligos_invoces.value = model['formul_oligos_invoces'].value
        self.fin_oligos_invoces.value = model['fin_oligos_invoces'].value
        self.total_price_oligos_invoces.value = model['total_price_oligos_invoces'].value
        self.ag_grid.options['rowData'] = model['ag_grid'].options['rowData']
        self.ag_grid.update()
        self.invoce_content_tab.options['rowData'] = model['invoce_content_tab'].options['rowData']
        self.invoce_content_tab.update()

