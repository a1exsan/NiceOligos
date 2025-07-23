from nicegui import ui
import pandas as pd

import xwell_plate_unit as XWells

class oligosynth_panel():
    def __init__(self, backend_model):
        self.backend_model = backend_model
        self.model = {}

        with ui.row().classes("w-full"):
            self.synth_name_label = ui.label('Synthesis name')
            self.synth_number_label = ui.label('Synthesis number')

        with (ui.grid(columns=3).classes("w-full").style("grid-template-columns: 1200px 200px 1300px")):
            with ui.grid(columns=3).classes("w-full").style("grid-template-columns: 200px 200px 200px"):
                self.on_clear_plate_button = ui.button("Clear plate", color="#FF1000").classes('w-[200px]')
            ui.label('Cell12')
            with ui.row():
                self.on_new_oligomap_button = ui.button("New Map", color="#E1A100").classes('w-[200px]')
                #options = ['AutoComplete', 'NiceGUI', 'Awesome']
                self.oligomap_name_input = ui.input(label='', placeholder='Map name')
                self.oligomap_syn_number_input = ui.input(label='', placeholder='Synth number')
                self.on_save_oligomap = ui.button('Save synth map', color="#00a100").classes('w-[200px]')
                self.on_update_oligomap = ui.button('Update map', color="orange").classes('w-[200px]')
                self.on_update_oligo_orders = ui.button('Update order', color="red").classes('w-[200px]')
                self.progressbar_1 = ui.spinner(size='md', color='#FFA000')
                self.progressbar_1.visible = False

            self.xwells_obj = XWells.XWell_plate(self, self.backend_model)

            with ui.column().classes("w-full"):
                self.number_of_copies_oligo = ui.input('', value=1)
                self.on_add_sel_oligo_to_plate = ui.button('Add selected', color="#e1a000").classes('w-[200px]')
                self.on_del_sel_oligo_to_plate = ui.button('Del selected', color="#f11000").classes('w-[200px]')
                ui.label('Select synthesis scale:')
                self.synth_scale_selector = ui.radio(['1 mg', '3 mg', '5 mg'], value='3 mg').props('inline')
                self.on_generate_oligomap_button = ui.button("Generate map", color="#00a100").classes('w-[200px]')
                with ui.dropdown_button('Set done operation', auto_close=True).classes('w-[200px]') as self.on_sel_done_btn:
                    for key in ['synth', 'cart', 'hplc', 'sed', 'paag', 'click', 'subl', 'Wasted', 'LCMS']:
                        ui.item(key, on_click=lambda n=key: self.done_event(n))
                self.on_sel_done_btn.props['color'] = 'red'
                with ui.dropdown_button('Set do operation', auto_close=True).classes('w-[200px]') as self.on_sel_do_btn:
                    for key in ['synth', 'cart', 'hplc', 'sed', 'paag', 'click', 'subl', 'Wasted', 'LCMS']:
                        ui.item(key, on_click=lambda n=key: self.do_event(n))
                self.on_sel_do_btn.props['color'] = 'orange'
                self.wells_layer_selector = ui.radio(
                    ['Base layer', 'Status layer', 'Purification layer', 'Support layer', 'Click layer'],
                    value='Base layer').props('inline')

            self.get_oligos_stack_grid()

        with ui.grid(columns=3).classes("w-full").style("grid-template-columns: 1200px 200px 1300px"):
            self.get_oligomaps_list_grid()
            with ui.column().classes("w-full"):
                self.on_show_actual_oligomaps = ui.button('Show actual maps', color="#E1A100").classes('w-[200px]')
                self.on_show_oligomaps = ui.button('Show all maps').classes('w-[200px]')
                self.on_load_oligomap = ui.button('Load map', color="#00a100").classes('w-[200px]')
            self.get_accord_tab()

        self.get_oligomap_grid()


    def done_event(self, e):
        print(e)

    def do_event(self, e):
        print(e)


    def get_model(self):
        self.model['xwells_obj'] = self.xwells_obj.get_copy()
        self.model['oligomap_ag_grid'] = self.oligomap_ag_grid
        self.model['oligomap_db_tab'] = self.oligomap_db_tab
        self.model['accord_tab'] = self.accord_tab
        return self.model


    def set_model(self, model):
        self.xwells_obj.set_copy(model['xwells_obj'])
        self.oligomap_ag_grid.options['rowData'] = model['oligomap_ag_grid'].options['rowData']
        self.oligomap_ag_grid.update()
        self.oligomap_db_tab.options['rowData'] = model['oligomap_db_tab'].options['rowData']
        self.oligomap_db_tab.update()
        self.oligomap_rowdata = model['oligomap_db_tab'].options['rowData']
        self.accord_tab.options['rowData'] = model['accord_tab'].options['rowData']
        self.accord_tab.update()


    def get_element_list(self, key):
        if key == 'button':
            return [
                ('on_new_oligomap_button', 'click'),
                ('on_add_sel_oligo_to_plate', 'click'),
                ('on_clear_plate_button', 'click'),
                ('on_generate_oligomap_button', 'click'),
                ('on_show_oligomaps', 'click'),
                ('on_show_actual_oligomaps', 'click'),
                ('on_load_oligomap', 'click'),
                ('on_del_sel_oligo_to_plate', 'click'),
                ('on_save_oligomap', 'click'),
                ('on_update_oligomap', 'click'),
                ('on_update_oligo_orders', 'click'),
            ]
        else:
            return []


    def __getitem__(self, item):
        if item == 'on_new_oligomap_button':
            return self.on_new_oligomap_button
        if item == 'on_add_sel_oligo_to_plate':
            return self.on_add_sel_oligo_to_plate
        if item == 'on_clear_plate_button':
            return self.on_clear_plate_button
        if item == 'on_generate_oligomap_button':
            return self.on_generate_oligomap_button
        if item == 'on_show_oligomaps':
            return self.on_show_oligomaps
        if item == 'on_show_actual_oligomaps':
            return self.on_show_actual_oligomaps
        if item == 'on_load_oligomap':
            return self.on_load_oligomap
        if item == 'on_del_sel_oligo_to_plate':
            return self.on_del_sel_oligo_to_plate
        if item == 'synth_scale_selector':
            return self.synth_scale_selector
        if item == 'on_save_oligomap':
            return self.on_save_oligomap
        if item == 'on_update_oligomap':
            return self.on_update_oligomap
        if item == 'on_sel_done_btn':
            return self.on_sel_done_btn
        if item == 'on_sel_do_btn':
            return self.on_sel_do_btn
        if item == 'on_update_oligo_orders':
            return self.on_update_oligo_orders
        if item == 'wells_layer_selector':
            return self.wells_layer_selector


    def get_oligomap_grid(self):
        map_tab = pd.DataFrame(
            {
                '#': [0],
                'map #': [0],
                'Order id': [0],
                'Position': [''],
                'Name': [''],
                'Sequence': [''],
                'Purif type': [''],
                'Support type': [''],
                'Date': [''],
                'Synt number': [''],
                'Scale, OE': [''],
                'CPG, mg': [''],
                'asm Sequence': [''],
                'Status': ['in queue'],
                'Dens, oe/ml': [0.],
                'Vol, ml': [0.3],
                'Purity, %': [50.],

                'Do LCMS': [True],
                'Done LCMS': [False],
                'Do synth': [True],
                'Done synth': [False],
                'Do cart': [True],
                'Done cart': [False],
                'Do hplc': [True],
                'Done hplc': [False],
                'Do paag': [True],
                'Done paag': [False],
                'Do click': [True],
                'Done click': [False],
                'Do sed': [True],
                'Done sed': [False],
                'Do subl': [True],
                'Done subl': [False],
                'DONE': [False],
                'Wasted': [False],
                'Send': [False],
            }
        )

        self.oligomap_rowdata = map_tab.to_dict('records')

        columnDefs = [
            {
                "field": "#",
                "checkboxSelection": True,
                "headerCheckboxSelection": True,
                "headerCheckboxSelectionFilteredOnly": True,
            },
            {"field": "map #"},
            {"field": "Order id"},
            {"field": "Position", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Name"},
            {"field": "Sequence", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Purif type", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Support type", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Date", 'filter': 'agTextColumnFilter'},
            {"field": "Synt number", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Scale, OE", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "CPG, mg", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "asm Sequence", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Status", 'editable': True, 'filter': 'agTextColumnFilter'},
            {"field": "Dens, oe/ml", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Vol, ml", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Purity, %", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},

            {"field": "Do LCMS", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done LCMS", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do synth", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done synth", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do cart", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done cart", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do hplc", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done hplc", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do paag", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done paag", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do sed", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done sed", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do click", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done click", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Do subl", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Done subl", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "DONE", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Wasted", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Send", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
        ]

        # with ui.column():
        self.oligomap_ag_grid = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': self.oligomap_rowdata,
                'rowSelection': 'multiple',
                "pagination": True,
                # "enableRangeSelection": True,
            }
            ,
            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
        self.oligomap_ag_grid.auto_size_columns = True
        self.oligomap_ag_grid.on("cellValueChanged", self.update_oligomap_cell_data)

    def update_oligomap_cell_data(self, e):
        self.oligomap_rowdata[e.args["rowIndex"]] = e.args["data"]

    def get_oligos_stack_grid(self):
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

        self.oligos_stack_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': invoce_content_df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                # "enterNavigatesVertically": True,
                # "enterNavigatesVerticallyAfterEdit": True,
                # "singleClickEdit": True
                ':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                                ' { background: "green" }',
            }
            ,
            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
        self.oligos_stack_tab.auto_size_columns = True

    def get_oligomaps_list_grid(self):
        map_db_tab = pd.DataFrame(
            {
                '#': [0],
                'Map name': [''],
                'Synth number': [''],
                'Date': [''],
                'in progress': [0],
                'Wasted': [0],
                'Finished': [0],
                'Total': [0],
                'avg yield, %': [0],
                'avg stage yield, %': [0]
            }
        )

        columnDefs = [
            {"field": "#"},
            {"field": "Map name"},
            {"field": "Synth number"},
            {"field": "Date"},
            {"field": "in progress"},
            {"field": "Wasted"},
            {"field": "Finished"},
            {"field": "Total"},
            {"field": "avg yield, %"},
            {"field": "avg stage yield, %"},
        ]

        self.oligomap_db_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': map_db_tab.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                # "enterNavigatesVertically": True,
                # "enterNavigatesVerticallyAfterEdit": True,
                # "singleClickEdit": True
                #':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                #                ' { background: "green" }',
            }
            ,
            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
        self.oligomap_db_tab.auto_size_columns = True

    def get_accord_tab(self):
        accord_tab = pd.DataFrame(
            {
                'Modification': ['A', 'C', 'G', 'T', '+A', '+C', '+G', '+T', '6FAM', 'Alk', 'R6G', 'HEX',
                                 'DEBL', 'ACTIV', 'CAPA', 'CAPB', 'OXID', 'R2', 'W1', 'W2'],
                'asm2000 position': ['A', 'C', 'G', 'T', '5', '6', '7', '8', '9', '[10]', '[11]', '[12]',
                                     'DEBL', 'ACTIV', 'CAPA', 'CAPB', 'OXID', 'R2', 'W1', 'W2'],
                'Conc, g/ml': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                               1, 1, 1, 1, 1, 1, 1, 1],
                'ul on step': [40., 40., 40., 40., 75., 75., 75., 75., 75., 75., 75., 75.,
                               240., 54., 45., 45., 110., 91., 110., 650.],
                'Amount 5mg, g': [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0.],
                'Amount 5mg, ml': [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                   0., 0., 0., 0., 0., 0., 0., 0.]
            }
        )

        columnDefs = [
            {"field": "Modification", 'editable': True},
            {"field": "asm2000 position", 'editable': True},
            {"field": "Conc, g/ml", 'editable': True},
            {"field": "ul on step", 'editable': True},
            {"field": "Amount 5mg, g", 'editable': True},
            # {"field": "ul on step, 10mg", 'editable': True},
            # {"field": "Amount 5mg, g"},
            {"field": "Amount 5mg, ml"},
            # {"field": "Amount 10mg, g"},
            # {"field": "Amount 10mg, ml"}
        ]

        self.accord_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': accord_tab.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True,
                # "enterNavigatesVertically": True,
                # "enterNavigatesVerticallyAfterEdit": True,
                # "singleClickEdit": True
                #':getRowStyle': '(params) => params.data.sufficiency < 0 ? { background: "red" } :'
                #                ' { background: "green" }',
            }
            ,
            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
        self.accord_tab.auto_size_columns = True

