from oligoMass import molmassOligo as mmo
from OligoMap_utils import api_db_interface
import pandas as pd
from nicegui import app, ui
from OligoMap_utils import oligomaps_search
import requests
from collections import Counter
from io import BytesIO

import xwell_plate_unit as XWells


class stat_dialog():
    def __init__(self, text1, text2):
        self.text1, self.text2 = text1, text2
        with ui.dialog() as self.dialog:
            with ui.card():
                ui.label(text=self.text1).style('font-size: 20px;')
                ui.label(text=self.text2).style('font-size: 20px;')
                with ui.row():
                    ui.button('Закрыть',  on_click=self.dialog.close).classes('w-[200px]')


class confirm_dialog():
    def __init__(self, text):
        self.text = text
        with ui.dialog() as self.dialog:
            with ui.card():
                ui.label(text=self.text).style('font-size: 20px;')
                with ui.row():
                    ui.button('Да', color='orange', on_click=self.do_confirm).classes('w-[200px]')
                    ui.button('Нет',  on_click=self.dialog.close).classes('w-[200px]')


    def on_confirm(self):
        pass

    def do_confirm(self):
        self.on_confirm()
        self.dialog.close()


class set_param_dialog():
    def __init__(self, rowdata, selrowdata):
        self.selrowdata = selrowdata
        self.rowdata = rowdata

        support_type = ['biocomma_500', 'biocomma_1000', 'biocomma_2000', 'biocomma_3000',
                        'bhq1_1000_hg', 'bhq2_1000_hg', 'bhq1_1000', 'bhq2_1000',
                        'bhq3_500', 'PO3_500']
        cpg_amount = ['1 mg', '2 mg', '3 mg', '4 mg', '5 mg', '6 mg', '7 mg', '8 mg']

        with ui.dialog() as self.dialog:
            with ui.card():
                ui.label('Установите параметры').style('font-size: 20px;')
                ui.label(f'Выбрано строк: {len(self.selrowdata)} ').style('font-size: 20px;')

                self.set_support = ui.select(options=support_type, with_input=True, label='Support type',
                            on_change=self.on_set_support_type_event).classes('w-[400px]').style('font-size: 20px;')

                self.set_cpg = ui.select(options=cpg_amount, with_input=True, label='Support amount',
                                            on_change=self.on_set_cpg_amount_event).classes('w-[400px]').style(
                    'font-size: 20px;')

                with ui.row():
                    self.volume = ui.input(label='volume').classes('w-[200px]').style('font-size: 20px;')
                    self.volume_btn = ui.button('Set volume', color='green', on_click=self.on_set_volume_event)

                with ui.row():
                    self.syn_number = ui.input(label='synth number').classes('w-[200px]').style('font-size: 20px;')
                    self.syn_number_btn = ui.button('Set synth number', color='green',
                                                    on_click=self.on_set_syn_number_event)


                with ui.row():
                    ui.button('Подтвердить', on_click=self.on_compleet_settings)
                    ui.button('Отмена', on_click=self.dialog.close)


    def set_param(self, param, value):
            sel_index, result = [], []
            if len(self.selrowdata) > 0:
                sel_index = list(pd.DataFrame(self.selrowdata)['#'])
            for row in self.rowdata:
                d = row.copy()
                if row['#'] in sel_index:
                    d[param] = value
                result.append(d)
            return result

    def on_set_volume_event(self):
        self.rowdata = self.set_param('Vol, ml', self.volume.value)

    def on_set_syn_number_event(self):
        self.rowdata = self.set_param('Synt number', self.syn_number.value)

    def on_set_support_type_event(self):
        self.rowdata = self.set_param('Support type', self.set_support.value)

    def on_set_cpg_amount_event(self):
        self.rowdata = self.set_param('CPG, mg', self.set_cpg.value)

    def on_compleet_settings(self):
        self.on_send_data(self.rowdata)
        self.dialog.close()

    def on_send_data(self, data):
        pass





class oligosynth_panel_page_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        self.pincode = app.storage.user.get('pincode')

        self.strftime_format = "%Y-%m-%d"
        self.db_name = 'scheduler_oligolab_2.db'

        with ui.row().classes("w-full"):
            self.synth_name_label = ui.label('Synthesis name')
            self.synth_number_label = ui.label('Synthesis number')

        with (ui.grid(columns=3).classes("w-full").style("grid-template-columns: 1200px 200px 1300px")):
            with ui.row().classes('gap-20'):
                self.on_clear_plate_button = ui.button("Clear plate", on_click=self.on_clear_plate_button_event,
                                                       color="#FF1000").classes('w-[200px]')
                self.on_selprint_excel_button = ui.button("Print label",
                                                on_click=self.on_selprint_excel_button_event).classes('w-[200px]')
                self.set_param_btn = ui.button("Set params", color='green',
                                                          on_click=self.on_set_param_btn_event).classes(
                    'w-[200px]')
            ui.label('Cell12')
            with ui.row():
                self.on_new_oligomap_button = ui.button("New Map", on_click=self.on_new_oligomap_button_event,
                                                        color="#E1A100").classes('w-[200px]')
                #options = ['AutoComplete', 'NiceGUI', 'Awesome']
                self.oligomap_name_input = ui.input(label='', placeholder='Map name')
                self.oligomap_syn_number_input = ui.input(label='', placeholder='Synth number')
                self.on_save_oligomap = ui.button('Save synth map',
                                on_click=self.on_save_oligomap_event, color="#00a100").classes('w-[200px]')
                self.on_update_oligomap = ui.button('Update map', on_click=self.on_update_oligomap_event,
                                                    color="orange").classes('w-[200px]')
                self.on_update_oligo_orders = ui.button('Update order', on_click=self.on_update_oligo_orders_event,
                                                        color="red").classes('w-[200px]')
                self.progressbar_1 = ui.spinner(size='md', color='#FFA000')
                self.progressbar_1.visible = False

            self.xwells_obj = XWells.XWell_plate()

            with ui.column().classes("w-full"):
                self.number_of_copies_oligo = ui.input('', value=1)
                self.on_add_sel_oligo_to_plate = ui.button('Add selected',
                                                           on_click=self.on_add_sel_oligo_to_plate_event,
                                                           color="#e1a000").classes('w-[200px]')
                self.on_del_sel_oligo_to_plate = ui.button('Del selected',
                                                            on_click=self.on_del_sel_oligo_to_plate_event,
                                                color="#f11000").classes('w-[200px]')
                self.on_replace_sel_oligo_to_plate = ui.button('Replace to selected',
                                                               on_click=self.on_replace_sel_oligo_to_plate_event,
                                                               color="#e1a000").classes('w-[200px]')
                ui.label('Select synthesis scale:')
                self.synth_scale_selector = ui.radio(['1 mg', '3 mg', '5 mg'], value='3 mg',
                                                     on_change=self.synth_scale_selector_event).props('inline')
                self.on_generate_oligomap_button = ui.button("Generate map",
                                                             on_click=self.on_generate_oligomap_button_event,
                                                             color="#00a100").classes('w-[200px]')
                with ui.dropdown_button('Set done operation', auto_close=True).classes('w-[200px]') as self.on_sel_done_btn:
                    for key in ['synth', 'cart', 'hplc', 'sed', 'paag', 'click', 'subl', 'Wasted', 'LCMS']:
                        ui.item(key, on_click=lambda n=key: self.on_sel_done_btn_event(n))
                self.on_sel_done_btn.props['color'] = 'red'
                with ui.dropdown_button('Set do operation', auto_close=True).classes('w-[200px]') as self.on_sel_do_btn:
                    for key in ['synth', 'cart', 'hplc', 'sed', 'paag', 'click', 'subl', 'Wasted', 'LCMS']:
                        ui.item(key, on_click=lambda n=key: self.on_sel_do_btn_event(n))
                self.on_sel_do_btn.props['color'] = 'orange'
                self.wells_layer_selector = ui.radio(
                    ['Base layer', 'Status layer', 'Purification layer', 'Support layer', 'Click layer',
                            'Order layer'],
                    on_change=self.wells_layer_selector_event,
                    value='Base layer').props('inline')

            self.get_oligos_stack_grid()

        with ui.grid(columns=3).classes("w-full").style("grid-template-columns: 1200px 200px 1300px"):
            self.get_oligomaps_list_grid()
            with ui.column().classes("w-full"):
                self.on_show_actual_oligomaps = ui.button('Show actual maps',
                                                          on_click=self.on_show_actual_oligomaps_event,
                                                            color="#E1A100").classes('w-[200px]')
                self.on_show_oligomaps = ui.button('Show all maps',
                                                   on_click=self.on_show_oligomaps_event).classes('w-[200px]')
                self.on_load_oligomap = ui.button('Load map', on_click=self.on_load_oligomap_event,
                                                  color="#00a100").classes('w-[200px]')
                self.delete_oligomap_btn = ui.button('Delete map', on_click=self.on_delete_oligomap_event,
                                                  color="red").classes('w-[200px]')
                self.delete_oligomap_btn = ui.button('Show stat', on_click=self.on_show_stat_event,
                                                     ).classes('w-[200px]')
            self.get_accord_tab()

        self.get_oligomap_grid()
        self.init_data()

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

    def update_accordtab_cell_data(self, e):
        self.accordtab_rowdata[e.args["rowIndex"]] = e.args["data"]
        self.accord_tab.options['rowData'] = self.accordtab_rowdata
        self.accord_tab.update()

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

    def get_init_accordtab_rowdata(self):
        return pd.DataFrame(
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
        ).to_dict('records')

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

        self.accordtab_rowdata = accord_tab.to_dict('records')

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
        self.accord_tab.on("cellValueChanged", self.update_accordtab_cell_data)


    def get_actual_stat_maps(self):
        ret = {'wasted %': 0., 'total wells': 0}
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        tab = oserch.get_oligomaps()
        orders = pd.DataFrame(self.get_orders_by_status('finished'))
        if len(tab) > 0:
            df = pd.DataFrame(tab)
            ret['wasted %'] = f"Wasted %: {round(df['Wasted'].sum() * 100/ df['Total'].sum(), 0)}"
            ret['total wells'] = f"Total wells: {df['Total'].sum()} / Total oligos: {orders.shape[0]}"
        return ret

    def on_show_actual_oligomaps_event(self):
        self.pincode = app.storage.user.get('pincode')
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData = oserch.get_actual_maps()
        self.oligomap_db_tab.options['rowData'] = rowData
        self.oligomap_db_tab.update()
        app.storage.user['oligomap_db_tab_rowdata'] = self.oligomap_db_tab.options['rowData']


    async def on_load_oligomap_event(self):
        self.pincode = app.storage.user.get('pincode')

        selRows = await self.oligomap_db_tab.get_selected_rows()
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData, accordData, synt_name, synt_number = oserch.load_oligomap(selRows)

        self.oligomap_ag_grid.options['rowData'] = rowData
        self.oligomap_ag_grid.update()
        self.oligomap_rowdata = rowData

        self.accord_tab.options['rowData'] = accordData
        self.accord_tab.update()

        self.synth_name_label.text = str(synt_name)
        self.synth_number_label.text = str(synt_number)

        self.xwells_obj.load_selrows(rowData)

        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']
        app.storage.user['accord_tab_rowdata'] = self.accord_tab.options['rowData']
        app.storage.user['synth_name_label'] = self.synth_name_label.text
        app.storage.user['synth_number_label'] = self.synth_number_label.text
        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()

    def delete_oligomap_event(self):
        selRows = app.storage.user.get('delete_oligomap_rows')
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        oserch.delete_map_from_base(selRows)
        self.on_show_actual_oligomaps_event()

    async def on_delete_oligomap_event(self):
        selRows = await self.oligomap_db_tab.get_selected_rows()
        app.storage.user['delete_oligomap_rows'] = selRows
        confirm = confirm_dialog(f'Удалить карту???')
        confirm.on_confirm = self.delete_oligomap_event
        confirm.dialog.open()


    def on_show_stat_event(self):
        stat_data = self.get_actual_stat_maps()
        #ret['wasted %'] = f"Wasted %: {round(df['Wasted'].sum() * 100 / df['Total'].sum(), 0)}"
        #ret['total wells']
        stat = stat_dialog(stat_data['wasted %'], stat_data['total wells'])
        stat.dialog.open()


    async def on_add_sel_oligo_to_plate_event(self):
        self.pincode = app.storage.user.get('pincode')
        selRows = await self.oligos_stack_tab.get_selected_rows()
        self.xwells_obj.add_selrows(selRows, int(self.number_of_copies_oligo.value))

        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def on_new_oligomap_button_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.oligos_stack_tab.options['rowData'] = []
        self.oligos_stack_tab.update()
        self.oligomap_rowdata = []
        self.accord_tab.options['rowData'] = self.get_init_accordtab_rowdata()
        self.accordtab_rowdata = self.get_init_accordtab_rowdata()
        self.accord_tab.update()
        self.oligomap_ag_grid.options['rowData'] = []
        self.oligomap_ag_grid.update()
        self.synth_name_label.text = ''
        self.synth_number_label.text = ''
        self.xwells_obj.clear_wells()

        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']
        app.storage.user['accord_tab_rowdata'] = self.accord_tab.options['rowData']
        app.storage.user['synth_name_label'] = self.synth_name_label.text
        app.storage.user['synth_number_label'] = self.synth_number_label.text
        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()
        app.storage.user['oligos_list_synth'] = []


    def on_clear_plate_button_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.xwells_obj.clear_wells()

        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def on_del_sel_oligo_to_plate_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.xwells_obj.clear_selected_wells()

        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def get_orders_by_status(self, status):

        def get_in_progress(find_list = ['synthesis', 'purification', 'formulation']):
            out = []
            for st in find_list:
                url = f'{self.api_db_url}/get_orders_by_status/{self.db_name}/{st}'
                ret = requests.get(url, headers=self.headers())
                out.extend(ret.json())
            return out
        out = []
        if status == 'in progress':
            out = get_in_progress(find_list = ['synthesis', 'purification', 'formulation'])
        elif status == 'total data':
            out = get_in_progress(find_list = ['in queue', 'synthesis', 'purification', 'formulation', 'finished'])
        else:
            url = f'{self.api_db_url}/get_orders_by_status/{self.db_name}/{status}'
            ret = requests.get(url, headers=self.headers())
            out.extend(ret.json())
        return out


    def wells_layer_selector_event(self, e):
        self.pincode = app.storage.user.get('pincode')

        if e.value == 'Order layer':
            orders = self.get_orders_by_status('total data')
            self.xwells_obj.set_invoce_data(orders)

        self.xwells_obj.layer_selector = e.value
        self.xwells_obj.draw_layers(e.value)

        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def on_replace_sel_oligo_to_plate_event(self):
        self.pincode = app.storage.user.get('pincode')
        self.xwells_obj.replace_to_selected_wells()

        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def save_passport(self, filename, data):
        buffer = BytesIO()
        data.to_excel(buffer, index=False)
        buffer.seek(0)
        ui.download(buffer.read(), filename=f'{filename}.xlsx')

    def on_selprint_excel_button_event(self):
        self.pincode = app.storage.user.get('pincode')

        pass_df = self.xwells_obj.get_selected_labels()
        self.save_passport('plate_label', pass_df)


    def return_scale_accord_tab(self, rowdata, scale):
        self.scale_dict = {'1 mg': 34., '3 mg': 40., '5 mg': 54.}
        out = []
        if self.check_pincode():
            for row in rowdata:
                d = row.copy()
                if row['Modification'] in 'A C G T'.split(' '):
                    #d['ul on step, 5mg'] = self.scale_dict[scale]
                    d['ul on step'] = self.scale_dict[scale]
                out.append(d)
            return out
        return rowdata


    def synth_scale_selector_event(self, e):
        self.pincode = app.storage.user.get('pincode')

        self.accord_tab.options['rowData'] = self.return_scale_accord_tab(
            self.accord_tab.options['rowData'], e.value)
        self.accord_tab.update()

        app.storage.user['accord_tab_rowdata'] = self.accord_tab.options['rowData']
        app.storage.user['scale_selector_value'] = e.value


    def on_save_oligomap_event(self):
        self.pincode = app.storage.user.get('pincode')
        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode
        map_name = self.oligomap_name_input.value
        map_synt_num = self.oligomap_syn_number_input.value
        rowData = self.oligomap_ag_grid.options['rowData']
        accord_rowData = self.accord_tab.options['rowData']
        omap.insert_map_to_base(map_name, map_synt_num, rowData, accord_rowData)
        self.on_show_actual_oligomaps_event()


    def on_update_oligomap_event(self):
        self.pincode = app.storage.user.get('pincode')

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        self.oligomap_ag_grid.options['rowData'] = self.oligomap_rowdata
        self.oligomap_ag_grid.update()
        rowData = self.oligomap_ag_grid.options['rowData']
        accord_rowData = self.accord_tab.options['rowData']
        rowData = omap.update_oligomap_status(rowData, accord_rowData)
        self.oligomap_ag_grid.options['rowData'] = rowData
        self.oligomap_ag_grid.update()
        self.xwells_obj.load_selrows(rowData)

        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']
        app.storage.user['xwells_obj'] = self.xwells_obj.get_copy()


    def on_update_oligo_orders_event(self):
        self.pincode = app.storage.user.get('pincode')

        self.progressbar_1.visible = True
        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode
        pos_list = self.xwells_obj.get_selected_pos_list()
        rowData_df = pd.DataFrame(self.oligomap_ag_grid.options['rowData'])
        sel_rowData_df = rowData_df[rowData_df['Position'].isin(pos_list)]
        omap.update_order_status(self.oligomap_ag_grid.options['rowData'], sel_rowData_df.to_dict('records'))

        self.progressbar_1.visible = False


    def on_show_oligomaps_event(self):
        self.pincode = app.storage.user.get('pincode')

        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData = oserch.get_oligomaps()
        self.oligomap_db_tab.options['rowData'] = rowData
        self.oligomap_db_tab.update()

        app.storage.user['oligomap_db_tab_rowdata'] = self.oligomap_db_tab.options['rowData']


    def on_sel_done_btn_event(self, e):
        self.pincode = app.storage.user.get('pincode')

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        rowData_df = pd.DataFrame(self.oligomap_ag_grid.options['rowData'])

        selected_wells = self.xwells_obj.selected_wells
        pos_list = []
        for key, well in zip(selected_wells.keys(), selected_wells.values()):
            pos_list.append(f"{well.symb}{well.num}")
        conditions = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))&(rowData_df[f"Done {e}"] == False)
        conditions_1 = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))&(rowData_df[f"Done {e}"] == True)
        rowData_df.loc[conditions, f"Done {e}"] = True
        rowData_df.loc[conditions_1, f"Done {e}"] = False

        self.oligomap_rowdata = omap.set_omap_status(rowData_df.to_dict('records'))
        self.oligomap_ag_grid.options['rowData'] = self.oligomap_rowdata
        self.oligomap_ag_grid.update()

        self.on_update_oligomap.run_method('click')
        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']


    def on_sel_do_btn_event(self, e):
        self.pincode = app.storage.user.get('pincode')

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        rowData_df = pd.DataFrame(self.oligomap_ag_grid.options['rowData'])

        selected_wells = self.xwells_obj.selected_wells
        pos_list = []
        for key, well in zip(selected_wells.keys(), selected_wells.values()):
            pos_list.append(f"{well.symb}{well.num}")
        conditions = (rowData_df[f"Do {e}"] == False)&(rowData_df['Position'].isin(pos_list))
        conditions_1 = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))
        rowData_df.loc[conditions, f"Do {e}"] = True
        rowData_df.loc[conditions_1, f"Do {e}"] = False

        self.oligomap_rowdata = omap.set_omap_status(rowData_df.to_dict('records'))
        self.oligomap_ag_grid.options['rowData'] = self.oligomap_rowdata
        self.oligomap_ag_grid.update()

        self.on_update_oligomap.run_method('click')
        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']



    def change_alk(self, rowData):
        out = []
        for row in rowData:
            oligo = mmo.oligoNASequence(row['Sequence'])
            tab = oligo.getSeqTabDF()
            #print(oligo.sequence)
            r = row.copy()
            if ((str(tab['prefix'].loc[1]).find('FAM') == -1) and (str(tab['prefix'].loc[1]) != '') and
                    (str(tab['prefix'].loc[1]).find('Alk') == -1)):
                if str(tab['prefix'].loc[1]).find('Biotin') == -1:
                    #sseq = [f'{i}{j}']
                    #seq = f"[Alk]{''.join(list(tab['nt']))}{str(tab['suffix'].loc[tab.shape[0]])}"
                    sseq = row['Sequence']
                    seq = sseq.replace(sseq[sseq.find('[') + 1: sseq.find(']')], 'Alk')
                    pref = str(tab['prefix'].loc[1])
                    pref = pref.replace('[', '')
                    pref = pref.replace(']', '')
                    purif_type = row['Purif type'] + f"_{pref}"
                    r['Sequence'] = seq
                    r['Purif type'] = purif_type
            out.append(r)

        return out


    def seq_to_asm_seq(self, accordRowData, tabRowData):
        tab = pd.DataFrame(tabRowData)
        accord = pd.DataFrame(accordRowData)

        seq_list = []
        for seq in tab['Sequence']:
            oligo = mmo.oligoNASequence(seq)
            df = oligo.getSeqTabDF()
            out_seq = ''
            for mod, nt in zip(df['prefix'], df['nt']):
                if '[' in mod and ']' in mod:
                    m = mod.replace('[', '')
                    m = m.replace(']', '')
                    out_seq += accord[accord['Modification'] == m]['asm2000 position'].max()
                    out_seq += accord[accord['Modification'] == nt]['asm2000 position'].max()
                elif mod == '+' or mod == '*' or mod == 'r' or mod == '':
                    m = f'{mod}{nt}'
                    if nt in 'R M S H V N Y K W B D N'.split(' '):
                        out_seq += nt
                    else:
                        out_seq += accord[accord['Modification'] == m]['asm2000 position'].max()
                else:
                    out_seq += nt

            seq_list.append(out_seq)
        tab['asm Sequence'] = seq_list
        tab['#'] = [i for i in range(1, len(seq_list) + 1)]
        return tab.to_dict('records')


    def get_all_amidites_types(self, seq_list):
        self.amidites_count = Counter()
        for seq in seq_list:
            #print(seq)
            o = mmo.oligoNASequence(seq)
            tab = o.getSeqTabDF()
            for mod, nt in zip(tab['prefix'], tab['nt']):
                if '[' in mod and ']' in mod and mod not in list(self.amidites_count.keys()):
                    self.amidites_count[mod] = 1
                    if nt not in list(self.amidites_count.keys()):
                        self.amidites_count[nt] = 1
                    else:
                        self.amidites_count[nt] += 1
                elif '[' in mod and ']' in mod and mod in list(self.amidites_count.keys()):
                    self.amidites_count[mod] += 1
                    if nt not in list(self.amidites_count.keys()):
                        self.amidites_count[nt] = 1
                    else:
                        self.amidites_count[nt] += 1
                elif mod == '+' or mod == '*' or mod == 'r' or mod == '':
                    a = f'{mod}{nt}'
                    if a not in list(self.amidites_count.keys()):
                        self.amidites_count[a] = 1
                    else:
                        self.amidites_count[a] += 1


    def update_accord_tab(self, accordRowData, tabRowData):
        tab = pd.DataFrame(tabRowData)
        accord = pd.DataFrame(accordRowData)

        round_n = 3
        constant_vol = 0.10

        self.get_all_amidites_types(list(tab['Sequence']))

        for mod, count in zip(self.amidites_count.keys(), self.amidites_count.values()):
            #print(mod, count)
            if '[' in mod and ']' in mod:
                m = mod.replace('[', '')
                m = m.replace(']', '')
                r5 = accord[accord['Modification'] == m]['ul on step'].max()
                r10 = 1#accord[accord['Modification'] == m]['ul on step, 10mg'].max()
                conc = accord[accord['Modification'] == m]['Conc, g/ml'].max()

                vol5 = r5 * count / 1000
                vol10 = r10 * count / 1000

                vol5 += vol5 * constant_vol
                vol10 += vol10 * constant_vol

                if vol5 <= 1.5:
                    vol5 = 1.5
                if vol10 <= 1.5:
                    vol10 = 1.5

                accord.loc[accord['Modification'] == m, 'Amount 5mg, ml'] = round(vol5, round_n)
                #accord.loc[accord['Modification'] == m, 'Amount 10mg, ml'] = round(vol10, round_n)
                accord.loc[accord['Modification'] == m, 'Amount 5mg, g'] = round(conc * vol5, round_n)
                #accord.loc[accord['Modification'] == m, 'Amount 10mg, g'] = round(conc * vol10, round_n)
            else:

                r5 = accord[accord['Modification'] == mod]['ul on step'].max()
                r10 = 1#accord[accord['Modification'] == mod]['ul on step, 10mg'].max()
                conc = accord[accord['Modification'] == mod]['Conc, g/ml'].max()

                vol5 = r5 * count / 1000
                vol10 = r10 * count / 1000

                vol5 += vol5 * constant_vol
                vol10 += vol10 * constant_vol

                if vol5 <= 1.5:
                    vol5 = 1.5
                if vol10 <= 1.5:
                    vol10 = 1.5

                accord.loc[accord['Modification'] == mod, 'Amount 5mg, ml'] = round(vol5, round_n)
                #accord.loc[accord['Modification'] == mod, 'Amount 10mg, ml'] = round(vol10, round_n)
                accord.loc[accord['Modification'] == mod, 'Amount 5mg, g'] = round(conc * vol5, round_n)
                #accord.loc[accord['Modification'] == mod, 'Amount 10mg, g'] = round(conc * vol10, round_n)

        total_count = sum(self.amidites_count.values())
        for mod, data in zip(accord['Modification'], accord['ul on step']):
            if mod in ['DEBL', 'ACTIV', 'CAPA', 'CAPB', 'OXID', 'R2', 'W1', 'W2']:
                accord.loc[accord['Modification'] == mod, 'Amount 5mg, ml'] = total_count * data / 1000

        return accord.to_dict('records')


    def on_generate_oligomap_button_event(self):
        self.pincode = app.storage.user.get('pincode')

        rowData = self.xwells_obj.get_rowData()
        rowData = self.change_alk(rowData)

        accord_rowdata = self.accord_tab.options['rowData']
        rowData = self.seq_to_asm_seq(accord_rowdata, rowData)
        accord_rowdata = self.update_accord_tab(accord_rowdata, rowData)

        self.accord_tab.options['rowData'] = accord_rowdata
        self.accord_tab.update()

        self.oligomap_ag_grid.options['rowData'] = rowData
        self.oligomap_ag_grid.update()

        self.oligomap_rowdata = rowData

        app.storage.user['accord_tab_rowdata'] = self.accord_tab.options['rowData']
        app.storage.user['oligomap_ag_grid_rowdata'] = self.oligomap_ag_grid.options['rowData']

    def on_set_param_btn_event(self):
        param_dialog = set_param_dialog(self.oligomap_ag_grid.options['rowData'], self.xwells_obj.get_selected_rows())
        param_dialog.on_send_data = self.on_set_param_dialog_rowdata
        param_dialog.dialog.open()

    def on_set_param_dialog_rowdata(self, data):
        self.oligomap_ag_grid.options['rowData'] = data
        self.oligomap_ag_grid.update()
        self.oligomap_rowdata = data
        self.on_update_oligomap.run_method('click')



    def init_data(self):
        if 'oligomap_db_tab_rowdata' in list(app.storage.user.keys()):
            self.oligomap_db_tab.options['rowData'] = app.storage.user.get('oligomap_db_tab_rowdata')
            self.oligomap_db_tab.update()

        if 'oligomap_ag_grid_rowdata' in list(app.storage.user.keys()):
            self.oligomap_ag_grid.options['rowData'] = app.storage.user.get('oligomap_ag_grid_rowdata')
            self.oligomap_ag_grid.update()

        if 'xwells_obj' in list(app.storage.user.keys()):
            self.xwells_obj.set_copy(app.storage.user.get('xwells_obj'))
            self.wells_layer_selector.value = self.xwells_obj.layer_selector

        if 'accord_tab_rowdata' in list(app.storage.user.keys()):
            self.accord_tab.options['rowData'] = app.storage.user.get('accord_tab_rowdata')
            self.accord_tab.update()

        if 'synth_name_label' in list(app.storage.user.keys()):
            self.synth_name_label.text = app.storage.user.get('synth_name_label')

        if 'synth_number_label' in list(app.storage.user.keys()):
            self.synth_number_label.text = app.storage.user.get('synth_number_label')

        if 'scale_selector_value' in list(app.storage.user.keys()):
            self.synth_scale_selector.value = app.storage.user.get('scale_selector_value')

        if 'oligos_list_synth' in list(app.storage.user.keys()):
            self.oligos_stack_tab.options['rowData'] = app.storage.user.get('oligos_list_synth')
            self.oligos_stack_tab.update()

