from oligoMass import molmassOligo as mmo
from OligoMap_utils import api_db_interface
import pandas as pd
from datetime import datetime
from nicegui import app
from OligoMap_utils import oligomaps_search

class Oligomap_backend(api_db_interface):
    def __init__(self, api_IP, db_port, stack):
        super().__init__(api_IP, db_port)

        self.oligomap_stack = stack
        self.strftime_format = "%Y-%m-%d"


    def init_frontend(self, front, ip, pincode):
        #print(pincode, self.client)
        if ip not in list(self.client.keys()):
            self.client[ip] = pincode
            self.client_frontend[ip] = front.get_model()
        self.client[ip] = pincode
        self.frontend = front
        self.frontend.set_model(self.client_frontend[ip])

        if ip in list(self.oligomap_stack.input_selected_rows.keys()):
            if len(self.oligomap_stack.input_selected_rows[ip]) > 0:
                self.frontend.oligos_stack_tab.options['rowData'] = self.oligomap_stack.input_selected_rows[ip][:96]
                self.frontend.oligos_stack_tab.update()

    def on_new_oligomap_button(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        #print(ip, len(self.oligomap_stack.input_selected_rows[ip]))
        if ip in list(self.oligomap_stack.input_selected_rows.keys()):
            self.oligomap_stack.input_selected_rows[ip] = []
            self.frontend.oligos_stack_tab.options['rowData'] = []
            self.frontend.oligos_stack_tab.update()


    async def on_add_sel_oligo_to_plate(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        #print(int(self.frontend.number_of_copies_oligo.value))
        selRows = await self.frontend.oligos_stack_tab.get_selected_rows()
        self.frontend.xwells_obj.add_selrows(selRows, int(self.frontend.number_of_copies_oligo.value))

        self.client_frontend[ip] = self.frontend.get_model()

    def on_clear_plate_button(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        self.frontend.xwells_obj.clear_wells()

        self.client_frontend[ip] = self.frontend.get_model()


    def on_generate_oligomap_button(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        rowData = self.frontend.xwells_obj.get_rowData()
        self.frontend.oligomap_ag_grid.options['rowData'] = rowData
        self.frontend.oligomap_ag_grid.update()

        self.client_frontend[ip] = self.frontend.get_model()


    def on_show_oligomaps(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData = oserch.get_oligomaps()
        self.frontend.oligomap_db_tab.options['rowData'] = rowData
        self.frontend.oligomap_db_tab.update()

        self.client_frontend[ip] = self.frontend.get_model()


    def on_show_actual_oligomaps(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData = oserch.get_actual_maps()
        self.frontend.oligomap_db_tab.options['rowData'] = rowData
        self.frontend.oligomap_db_tab.update()

        self.client_frontend[ip] = self.frontend.get_model()


    async def on_load_oligomap(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        selRows = await self.frontend.oligomap_db_tab.get_selected_rows()
        oserch = oligomaps_search(self.db_IP, self.db_port)
        oserch.pincode = self.pincode
        rowData, accordData, synt_name, synt_number = oserch.load_oligomap(selRows)

        self.frontend.oligomap_ag_grid.options['rowData'] = rowData
        self.frontend.oligomap_ag_grid.update()

        self.frontend.accord_tab.options['rowData'] = accordData
        self.frontend.accord_tab.update()

        self.frontend.synth_name_label.text = str(synt_name)
        self.frontend.synth_number_label.text = str(synt_number)

        self.frontend.xwells_obj.load_selrows(rowData)

        self.client_frontend[ip] = self.frontend.get_model()


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