from oligoMass import molmassOligo as mmo
from OligoMap_utils import api_db_interface
import pandas as pd
from nicegui import app, ui
from OligoMap_utils import oligomaps_search
import requests
from collections import Counter
from io import BytesIO

class Oligomap_backend(api_db_interface):
    def __init__(self, api_IP, db_port, stack):
        super().__init__(api_IP, db_port)

        self.oligomap_stack = stack
        self.strftime_format = "%Y-%m-%d"
        self.db_name = 'scheduler_oligolab_2.db'


    def init_frontend(self, front, ip, pincode):
        #print(pincode, self.client)
        if ip not in list(self.client.keys()):
            self.client[ip] = pincode
            self.client_frontend[ip] = front.get_model()

        #print(ip, pincode)
        self.client[ip] = pincode
        self.frontend = front
        self.frontend.set_model(self.client_frontend[ip])
        #print(self.client_frontend[ip])

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
        self.frontend.oligomap_rowdata = []
        self.frontend.accord_tab.options['rowData'] = self.frontend.get_init_accordtab_rowdata()
        self.frontend.accordtab_rowdata = self.frontend.get_init_accordtab_rowdata()
        self.frontend.accord_tab.update()
        self.frontend.oligomap_ag_grid.options['rowData'] = []
        self.frontend.oligomap_ag_grid.update()
        self.frontend.synth_name_label.text = ''
        self.frontend.synth_number_label.text = ''

        self.frontend.set_model()


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
        rowData = self.change_alk(rowData)

        accord_rowdata = self.frontend.accord_tab.options['rowData']
        rowData = self.seq_to_asm_seq(accord_rowdata, rowData)
        accord_rowdata = self.update_accord_tab(accord_rowdata, rowData)

        self.frontend.accord_tab.options['rowData'] = accord_rowdata
        self.frontend.accord_tab.update()

        self.frontend.oligomap_ag_grid.options['rowData'] = rowData
        self.frontend.oligomap_ag_grid.update()

        self.frontend.oligomap_rowdata = rowData

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
        self.frontend.oligomap_rowdata = rowData

        self.frontend.accord_tab.options['rowData'] = accordData
        self.frontend.accord_tab.update()

        self.frontend.synth_name_label.text = str(synt_name)
        self.frontend.synth_number_label.text = str(synt_number)

        self.frontend.xwells_obj.load_selrows(rowData)

        self.client_frontend[ip] = self.frontend.get_model()
        #print(self.client_frontend.keys(), ip)
        #print(self.client_frontend[ip]['xwells_obj'])


    def on_del_sel_oligo_to_plate(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.frontend.xwells_obj.clear_selected_wells()

        self.client_frontend[ip] = self.frontend.get_model()

    def on_replace_sel_oligo_to_plate(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.frontend.xwells_obj.replace_to_selected_wells()

        self.client_frontend[ip] = self.frontend.get_model()

    def synth_scale_selector(self, e):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.frontend.accord_tab.options['rowData'] = self.return_scale_accord_tab(
            self.frontend.accord_tab.options['rowData'], e.value)
        self.frontend.accord_tab.update()

        self.client_frontend[ip] = self.frontend.get_model()


    def wells_layer_selector(self, e):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]
        #ui.run_javascript('location.reload()')

        if e.value == 'Order layer':
            orders = self.get_orders_by_status('total data')
            self.frontend.xwells_obj.set_invoce_data(orders)

        self.frontend.xwells_obj.layer_selector = e.value
        self.frontend.xwells_obj.draw_layers(e.value)

        self.client_frontend[ip] = self.frontend.get_model()


    def on_save_oligomap(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        map_name = self.frontend.oligomap_name_input.value
        map_synt_num = self.frontend.oligomap_syn_number_input.value

        rowData = self.frontend.oligomap_ag_grid.options['rowData']
        accord_rowData = self.frontend.accord_tab.options['rowData']

        omap.insert_map_to_base(map_name, map_synt_num, rowData, accord_rowData)

        self.client_frontend[ip] = self.frontend.get_model()


    def on_update_oligomap(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        self.frontend.oligomap_ag_grid.options['rowData'] = self.frontend.oligomap_rowdata
        self.frontend.oligomap_ag_grid.update()
        rowData = self.frontend.oligomap_ag_grid.options['rowData']
        accord_rowData = self.frontend.accord_tab.options['rowData']
        rowData = omap.update_oligomap_status(rowData, accord_rowData)
        self.frontend.oligomap_ag_grid.options['rowData'] = rowData
        self.frontend.oligomap_ag_grid.update()
        self.frontend.xwells_obj.load_selrows(rowData)

        self.client_frontend[ip] = self.frontend.get_model()


    async def on_update_oligo_orders(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        self.frontend.progressbar_1.visible = True
        selrows = await self.frontend.oligomap_ag_grid.get_selected_rows()

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        pos_list = self.frontend.xwells_obj.get_selected_pos_list()
        rowData_df = pd.DataFrame(self.frontend.oligomap_ag_grid.options['rowData'])
        sel_rowData_df = rowData_df[rowData_df['Position'].isin(pos_list)]
        omap.update_order_status(self.frontend.oligomap_ag_grid.options['rowData'], sel_rowData_df.to_dict('records'))

        self.client_frontend[ip] = self.frontend.get_model()
        self.frontend.progressbar_1.visible = False


    def on_sel_done_btn(self, e):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        rowData_df = pd.DataFrame(self.frontend.oligomap_ag_grid.options['rowData'])

        selected_wells = self.frontend.xwells_obj.selected_wells
        pos_list = []
        for key, well in zip(selected_wells.keys(), selected_wells.values()):
            pos_list.append(f"{well.symb}{well.num}")
        conditions = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))&(rowData_df[f"Done {e}"] == False)
        conditions_1 = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))&(rowData_df[f"Done {e}"] == True)
        rowData_df.loc[conditions, f"Done {e}"] = True
        rowData_df.loc[conditions_1, f"Done {e}"] = False

        self.frontend.oligomap_rowdata = omap.set_omap_status(rowData_df.to_dict('records'))
        self.frontend.oligomap_ag_grid.options['rowData'] = self.frontend.oligomap_rowdata
        self.frontend.oligomap_ag_grid.update()

        self.frontend.on_update_oligomap.run_method('click')

        self.client_frontend[ip] = self.frontend.get_model()


    def on_sel_do_btn(self, e):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        omap = oligomaps_search(self.db_IP, self.db_port)
        omap.pincode = self.pincode

        rowData_df = pd.DataFrame(self.frontend.oligomap_ag_grid.options['rowData'])

        selected_wells = self.frontend.xwells_obj.selected_wells
        pos_list = []
        for key, well in zip(selected_wells.keys(), selected_wells.values()):
            pos_list.append(f"{well.symb}{well.num}")
        conditions = (rowData_df[f"Do {e}"] == False)&(rowData_df['Position'].isin(pos_list))
        conditions_1 = (rowData_df[f"Do {e}"] == True)&(rowData_df['Position'].isin(pos_list))
        rowData_df.loc[conditions, f"Do {e}"] = True
        rowData_df.loc[conditions_1, f"Do {e}"] = False

        self.frontend.oligomap_rowdata = omap.set_omap_status(rowData_df.to_dict('records'))
        self.frontend.oligomap_ag_grid.options['rowData'] = self.frontend.oligomap_rowdata
        self.frontend.oligomap_ag_grid.update()

        self.frontend.on_update_oligomap.run_method('click')

        self.client_frontend[ip] = self.frontend.get_model()


    def on_selprint_excel_button(self):
        ip = app.storage.user.get('client_ip')
        self.pincode = self.client[ip]

        pass_df = self.frontend.xwells_obj.get_selected_labels()
        self.save_passport('plate_label', pass_df)

        self.client_frontend[ip] = self.frontend.get_model()

    def check_pincode(self):
        url = f'{self.api_db_url}/get_all_invoces/{self.db_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.status_code == 200


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
        if item == 'on_selprint_excel_button':
            return self.on_selprint_excel_button
        if item == 'on_replace_sel_oligo_to_plate':
            return self.on_replace_sel_oligo_to_plate


    def save_passport(self, filename, data):
        buffer = BytesIO()
        data.to_excel(buffer, index=False)
        buffer.seek(0)
        ui.download(buffer.read(), filename=f'{filename}.xlsx')