
from nicegui import ui, app
import re
from chemicals_page import moleculeInfo
from Reactor import image_background
from OligoMap_utils import api_db_interface
import json
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D


class modification():
    def __init__(self, symbol, smiles=None, unicode=None, data_json=None):
        self.symbol = symbol
        self.smiles = smiles
        self.unicode = unicode
        self.data_json = data_json

    def get_mol_mass(self):
        if self.smiles is not None and self.smiles != '':
            mol = moleculeInfo(self.smiles)
            return mol.get_props()['Mol weight, Da']
        else:
            return 0.

    def to_dict(self):
        return self.__dict__

    @classmethod
    def from_dict(cls, dict_data):
        return cls(**dict_data)

    def draw_mod_svg(self, context, width=800, height=400):
        #context.content = ''
        if self.smiles != '':
            mol = moleculeInfo(self.smiles)
            svg = mol.draw_svg(width, height)
            context.content += f'<g transform="translate(50, {height + 100})">{svg}</g>'

    def draw_mod_svg_reactant(self, context, width=800, height=400):
        #context.content = ''
        if self.smiles != '':
            mol = moleculeInfo(self.smiles)
            svg = mol.draw_svg(width, height)
            context.content += f'<g transform="translate({width + 100}, {height + 100})">{svg}</g>'


class reaction_smarts():
    def __init__(self, name, smarts=None, data_json=None):
        self.name = name
        self.smarts = smarts
        self.data_json = data_json

    def to_dict(self):
        return self.__dict__

    @classmethod
    def from_dict(cls, dict_data):
        return cls(**dict_data)

    def draw_rnx_svg(self, context, width=800, height=400):
        #context.content = ''
        if self.smarts != '':
            reaction_obj = AllChem.ReactionFromSmarts(self.smarts)
            drawer = Draw.MolDraw2DSVG(width, height)
            rdMolDraw2D.SetDarkMode(drawer)
            drawer.DrawReaction(reaction_obj)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            context.content += f'<g transform="translate(50, 50)">{svg}</g>'



class modification_base(api_db_interface):
    def __init__(self, db_IP, db_port):
        super().__init__(db_IP, db_port)
        self.reaction_base = {}
        self.modification_base = {}
        self.rnx_sel_id = 0
        self.mod_sel_id = 0
        self.reactant = modification.from_dict({
            'symbol':'',
            'unicode':'',
            'smiles':'',
            'data_json':''
                                                })
        self.context_width = 1400
        self.context_height = 1200
        self.border_color = 'gray'


    def draw_context(self, context):
        context.content = ''
        context.content += (f'<rect x={0} y={0} '
                            f'width={self.context_width} height={self.context_height} '
                            f' rx=10 ry=10 fill="{None}" fill-opacity="{0.6}"'
                            f'stroke="{self.border_color}" stroke-width="4"/>')
        if self.rnx_sel_id > 0:
            self.reaction_base[self.rnx_sel_id].draw_rnx_svg(context, width=600, height=400)
        if self.mod_sel_id > 0:
            self.modification_base[self.mod_sel_id].draw_mod_svg(context, width=600, height=400)
            if self.reactant.smiles != '':
                self.reactant.draw_mod_svg_reactant(context, width=600, height=400)

    def get_all_data_from_base(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.mod_db_name}/{tab_name}'
        r = requests.get(url, headers=self.headers())
        if r.status_code == 200:
            return r.json()
        else:
            return []

    def get_reaction_rowdata(self):
        tab = self.get_all_data_from_base('rnx')
        out = []
        self.reaction_base = {}
        for row in tab:
            d = {}
            d['name'] = row[1]
            d['smarts'] = row[2]
            d['data_json'] = row[3]
            self.reaction_base[row[0]] = reaction_smarts.from_dict(d)
            d['id'] = row[0]
            out.append(d)
        return out

    def single_reaction(self):
        if self.reactant.smiles != '' and self.rnx_sel_id > 0:
            rxn = rdChemReactions.ReactionFromSmarts(self.reaction_base[self.rnx_sel_id].smarts)
            react = Chem.MolFromSmiles(self.reactant.smiles)
            products = rxn.RunReactants([react])
            if products != []:
                self.reactant.smiles = Chem.MolToSmiles(products[0][0])

    def get_modification_rowdata(self):
        tab = self.get_all_data_from_base('mod')
        out = []
        self.modification_base = {}
        for row in tab:
            d = {}
            d['symbol'] = row[1]
            d['unicode'] = row[2]
            d['smiles'] = row[3]
            d['data_json'] = row[4]
            self.modification_base[row[0]] = modification.from_dict(d)
            d['id'] = row[0]
            out.append(d)
        return out

    def insert_modification(self, data_object):
        tab = self.get_all_data_from_base('mod')
        insert = True
        for row in tab:
            if row[1] == data_object['symbol']:
                insert = False
                break
        if insert:
            url = f'{self.api_db_url}/insert_data/{self.mod_db_name}/mod'
            r = requests.post(url,
                          json=json.dumps([
                              data_object['symbol'],
                              data_object['unicode'],
                              data_object['smiles'],
                              data_object['data_json']
                          ]), headers=self.headers())
            if r.status_code == 200:
                ui.notify(f'Модификция добавлена в базу')
        else:
            ui.notify(f'Модификция есть в базе')

    def insert_reaction(self, data_object):
        tab = self.get_all_data_from_base('rnx')
        insert = True
        for row in tab:
            if row[1] == data_object['name']:
                insert = False
                break
        if insert:
            url = f'{self.api_db_url}/insert_data/{self.mod_db_name}/rnx'
            r = requests.post(url,
                          json=json.dumps([
                              data_object['name'],
                              data_object['smarts'],
                              data_object['data_json']
                          ]), headers=self.headers())
            if r.status_code == 200:
                ui.notify(f'Реакция добавлена в базу')
        else:
            ui.notify(f'Реакция есть в базе')

    def update_reaction(self, id, reaction):
        data = reaction.to_dict()
        url = f'{self.api_db_url}/update_data/{self.mod_db_name}/rnx/{id}'
        r = requests.put(url,
                         json=json.dumps({
                             'name_list': list(data.keys()),
                             'value_list': list(data.values())
                         }), headers=self.headers())
        if r.status_code == 200:
            ui.notify(f'Реакция успешно обновлена')

    def update_modification(self, id, mod):
        data = mod.to_dict()
        url = f'{self.api_db_url}/update_data/{self.mod_db_name}/mod/{id}'
        r = requests.put(url,
                         json=json.dumps({
                             'name_list': list(data.keys()),
                             'value_list': list(data.values())
                         }), headers=self.headers())
        if r.status_code == 200:
            ui.notify(f'Модификация успешно обновлена')



class single_nucleic_acid_chain():
    def __init__(self, seq):
        self.sequence = seq
        self.chain = self.parse(self.sequence)

    def parse(self, seq):
        seq = re.sub(r'\+(\w)', r'[+\1]', seq)
        pattern = r'(\[[^\]]*\])|([^\[\]]+)'
        matches = re.findall(pattern, seq)
        result = []
        for match in matches:
            bracket_text, normal_text = match
            if bracket_text:
                result.append(bracket_text)
            else:
                cleaned = normal_text.strip()
                if cleaned:
                    result.extend(list(cleaned.upper()))
        return result

class single_nucleic_acid_chain_assembler(single_nucleic_acid_chain):
    def __init__(self, seq, rnx_base, mod_base):
        super().__init__(seq)
        self.rnx_base = rnx_base
        base = {}
        for val in mod_base.values():
            base[val.symbol] = val
        self.mod_base = base
        self.structure = ''

    def do_auto_reactions(self, modification):
        if self.structure == '':
            self.structure = modification.smiles

        print(modification.symbol)
        print(self.structure)

        adduct = modification.smiles
        rnx_id_list = json.loads(modification.data_json)['rnx_id_list']

        for i, id in enumerate(rnx_id_list):
            print(i ,id)
            if i == 0:
                rxn = rdChemReactions.ReactionFromSmarts(self.rnx_base[id].smarts)
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react])
                print(products)
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])
            if i == 1:
                rxn = rdChemReactions.ReactionFromSmarts(self.rnx_base[id].smarts)
                react2 = Chem.MolFromSmiles(self.structure)
                react1 = Chem.MolFromSmiles(adduct)
                products = rxn.RunReactants([react1, react2])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])
                    print(self.structure)
            if i > 1:
                rxn = rdChemReactions.ReactionFromSmarts(self.rnx_base[id].smarts)
                if id == 4:
                    react = Chem.MolFromSmiles(self.structure)
                    ox = Chem.MolFromSmiles('O')
                    products = rxn.RunReactants([react, ox])
                else:
                    react = Chem.MolFromSmiles(self.structure)
                    products = rxn.RunReactants([react])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])
                    print(self.structure)




    def build(self):
        reverse_chain = self.chain[::-1]
        for token in reverse_chain:
            if token in self.mod_base:
                #print(self.mod_base[token].to_dict())
                self.do_auto_reactions(self.mod_base[token])
            else:
                print(f'{token} not in base')




def run():
    ss_dna = single_nucleic_acid_chain('[A]CGTA[FAM_dT[mod2]CGTIIUa+agggURSY+TTT[BHQ2]')
    print(ss_dna.chain)


    d = {'symbol': 'A',
         'smiles': 'CC(N([P@](O[C@@]1([C@@](O[C@](N2C3=NC=NC(/N=C(\O)/C4C=CC=CC=4)=C3N=C2)([H])C1)([H])COC(C1C=CC(OC)=CC=1)(C1C=CC(OC)=CC=1)C1C=CC=CC=1)[H])OCCC#N)C(C)C)C',
         'unicode': 'abracodabra'
         }
    mod = modification.from_dict(d)
    #mod = modification('A')
    #mod.smiles = 'CC(N([P@](O[C@@]1([C@@](O[C@](N2C3=NC=NC(/N=C(\O)/C4C=CC=CC=4)=C3N=C2)([H])C1)([H])COC(C1C=CC(OC)=CC=1)(C1C=CC(OC)=CC=1)C1C=CC=CC=1)[H])OCCC#N)C(C)C)C'
    print(mod.get_mol_mass())
    print(mod.to_dict())

def run_mods():
    IP = '127.0.0.1'
    port = '8012'
    pin = ''
    mod_base = modification_base(IP, port)
    mod_base.pincode = pin

    rnx_obj = {}
    rnx_obj['name'] = 'DETRIT DMT'
    rnx_obj['smarts'] = 'C-O-c1:c:c:c(-C(-[O;H0;D2;+0:1]-[C:2])(-c2:c:c:c:c:c:2)-c2:c:c:c(-O-C):c:c:2):c:c:1>>[C:2]-[OH;D1;+0:1]'
    rnx_obj['data_json'] = ''
    reaction = reaction_smarts.from_dict(rnx_obj)
    print(reaction.to_dict())
    mod_base.insert_reaction(reaction.to_dict())

    rnx_obj = {}
    rnx_obj['name'] = 'COUPLING main'
    rnx_obj[
        'smarts'] = '[#8:1]-[P@;H0;D3;+0:2](-[#8:3])-[N;H0;D3;+0:4](-[C:5])-[C:6].[C:7]-[OH;D1;+0:8]>>[#8:1]-[P;H0;D3;+0:2](-[#8:3])-[O;H0;D2;+0:8]-[C:7].[C:5]-[NH;D2;+0:4]-[C:6]'
    rnx_obj['data_json'] = ''
    #mod_base.insert_reaction(rnx_obj)

    rnx_obj = {}
    rnx_obj['name'] = 'DEBENZYL dA'
    rnx_obj[
        'smarts'] = 'O/C(=[N;H0;D2;+0:1]\[c:2](:[#7;a:3]):[c:4]:[#7;a:5])-c1:c:c:c:c:c:1>>[#7;a:5]:[c:4]:[c:2](:[#7;a:3])-[NH2;D1;+0:1]'
    rnx_obj['data_json'] = ''
    reaction = reaction_smarts.from_dict(rnx_obj)
    print(reaction.to_dict())
    mod_base.insert_reaction(reaction.to_dict())


class modification_page_model():
    def __init__(self):
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        self.obj_base = modification_base(IP, port)
        self.obj_base.pincode = app.storage.user.get('pincode')
        rowdata = self.obj_base.get_reaction_rowdata()
        mod_rowdata = self.obj_base.get_modification_rowdata()

        colDefs = [
            {"field": "id", 'editable': False},
            {"field": "name", 'editable': True},
            {"field": "smarts", 'editable': True},
            {"field": "data_json", 'editable': True}
        ]

        colDefs_modif = [
            {"field": "id", 'editable': False},
            {"field": "symbol", 'editable': True},
            {"field": "smiles", 'editable': True},
            {"field": "unicode", 'editable': True},
            {"field": "data_json", 'editable': True}
        ]

        with ui.row():
            with ui.column():
                ui.button('add rnx', on_click=self.on_add_reaction)
                ui.button('add mod', color='green', on_click=self.on_add_modification)
                ui.button('react', color='orange', on_click=self.on_react_modification)
                ui.button('edit rnx json', on_click=self.on_edit_rnx_json)
                ui.button('edit mod json', color='green', on_click=self.on_edit_mod_json)
            with ui.column():
                self.rnx_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': rowdata,
                        'rowSelection': 'multiple',
                        "pagination": True,
                        "enterNavigatesVertically": True,
                        "enterNavigatesVerticallyAfterEdit": True,
                        "singleClickEdit": True,
                        # "enableRangeSelection": True,
                    },
                    theme='alpine-dark').style('height: 600px; width: 1200px')  # alpine  material  quartz  balham
                self.rnx_grid.auto_size_columns = True
                self.rnx_grid.on("cellValueChanged", self.update_grid_cell_data)
                self.rnx_grid.on('rowSelected', self.on_select_rnx_row)

                self.mod_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs_modif,
                        'rowData': mod_rowdata,
                        'rowSelection': 'multiple',
                        "pagination": True,
                        "enterNavigatesVertically": True,
                        "enterNavigatesVerticallyAfterEdit": True,
                        "singleClickEdit": True,
                        # "enableRangeSelection": True,
                    },
                    theme='alpine-dark').style('height: 600px; width: 1200px')
                self.mod_grid.auto_size_columns = True
                self.mod_grid.on("cellValueChanged", self.update_grid_cell_data_modif)
                self.mod_grid.on('rowSelected', self.on_select_mod_row)

            with ui.column():
                self.img_width = 1350
                self.img_height = 1200
                self.bkg_color = '#010101'
                bkg = image_background(self.img_width, self.img_height, color=self.bkg_color)
                self.image = ui.interactive_image(f"data:image/png;base64,{bkg.background_base64}",
                                                  size=(self.img_width, self.img_height), on_mouse=self.mouse_handler,
                                                  events=['mousedown', 'mousemove', 'mouseup', 'click']
                                                  )
                self.rnx_context = self.image.add_layer()

        with ui.row():
            ui.button('Parce seq', on_click=self.on_parse_seq)
            self.seq_area = ui.textarea(label='Sequence').style('width: 400px')
            self.tokens_area = ui.textarea(label='Tokens').style('width: 400px')

            self.init_reaction_rowdata()

    def init_reaction_rowdata(self):
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        self.obj_base = modification_base(IP, port)
        self.obj_base.pincode = app.storage.user.get('pincode')
        self.rnx_grid.options['rowData'] = self.obj_base.get_reaction_rowdata()
        self.rnx_grid.update()
        self.mod_grid.options['rowData'] = self.obj_base.get_modification_rowdata()
        self.mod_grid.update()
        self.obj_base.context_width = self.img_width
        self.obj_base.context_height = self.img_height
        self.obj_base.draw_context(self.rnx_context)

    def update_grid_cell_data(self, e):
        d = {'name': e.args['data']['name'],
             'smarts': e.args['data']['smarts'],
             'data_json': e.args['data']['data_json']}
        reaction = reaction_smarts.from_dict(d)
        self.obj_base.update_reaction(e.args['data']['id'], reaction)
        self.init_reaction_rowdata()

    def update_grid_cell_data_modif(self, e):
        d = {'symbol': e.args['data']['symbol'],
             'unicode': e.args['data']['unicode'],
             'smiles': e.args['data']['smiles'],
             'data_json': e.args['data']['data_json']}
        mod = modification.from_dict(d)
        self.obj_base.update_modification(e.args['data']['id'], mod)
        self.init_reaction_rowdata()

    async def on_select_mod_row(self):
        selrows = await self.mod_grid.get_selected_rows()
        self.obj_base.mod_sel_id = selrows[0]['id']
        self.obj_base.draw_context(self.rnx_context)
        self.obj_base.reactant = modification.from_dict(self.obj_base.modification_base[self.obj_base.mod_sel_id].to_dict())


    def mouse_handler(self, e):
        pass

    async def on_select_rnx_row(self):
        selrows = await self.rnx_grid.get_selected_rows()
        self.obj_base.rnx_sel_id = selrows[0]['id']
        self.obj_base.draw_context(self.rnx_context)
        #self.obj_base.reaction_base[selrows[0]['id']].draw_rnx_svg(self.rnx_context, width=600, height=400)

    def on_add_reaction(self):
        rnx_obj = {}
        rnx_obj['name'] = 'default'
        rnx_obj[
            'smarts'] = ''
        rnx_obj['data_json'] = ''
        reaction = reaction_smarts.from_dict(rnx_obj)
        self.obj_base.insert_reaction(reaction.to_dict())
        self.init_reaction_rowdata()

    def on_add_modification(self):
        mod_obj = {}
        mod_obj['symbol'] = 'default'
        mod_obj[
            'unicode'] = ''
        mod_obj['smiles'] = ''
        mod_obj['data_json'] = ''
        mod = modification.from_dict(mod_obj)
        self.obj_base.insert_modification(mod.to_dict())
        self.init_reaction_rowdata()

    def on_react_modification(self):
        self.obj_base.single_reaction()
        self.obj_base.draw_context(self.rnx_context)

    async def on_edit_rnx_json(self):
        selrows = await self.rnx_grid.get_selected_rows()
        editor = edit_json_dialog(selrows[0]['data_json'], f'Edit reaction {selrows[0]["name"]}')
        editor.on_save = self.on_save_rnx_json
        editor.dialog.open()

    async def on_edit_mod_json(self):
        selrows = await self.mod_grid.get_selected_rows()
        editor = edit_json_dialog(selrows[0]['data_json'], f'Edit modification {selrows[0]["symbol"]}')
        editor.on_save = self.on_save_mod_json
        editor.dialog.open()

    def on_save_rnx_json(self, data):
        if self.obj_base.rnx_sel_id > 0:
            rnx_dict = self.obj_base.reaction_base[self.obj_base.rnx_sel_id].to_dict()
            d = {'name': rnx_dict['name'],
                 'smarts': rnx_dict['smarts'],
                 'data_json': data['json']}
            reaction = reaction_smarts.from_dict(d)
            self.obj_base.update_reaction(self.obj_base.rnx_sel_id, reaction)
            self.init_reaction_rowdata()

    def on_save_mod_json(self, data):
        if self.obj_base.mod_sel_id > 0:
            mod_dict = self.obj_base.modification_base[self.obj_base.mod_sel_id].to_dict()
            d = {'symbol': mod_dict['symbol'],
                 'unicode': mod_dict['unicode'],
                 'smiles': mod_dict['smiles'],
                 'data_json': data['json']}
            mod = modification.from_dict(d)
            self.obj_base.update_modification(self.obj_base.mod_sel_id, mod)
            self.init_reaction_rowdata()

    def on_parse_seq(self):
        sequence = self.seq_area.value
        oligo = single_nucleic_acid_chain(sequence)
        self.tokens_area.value = json.dumps(oligo.chain)

        oligo = single_nucleic_acid_chain_assembler(sequence,
                                                      self.obj_base.reaction_base,
                                                      self.obj_base.modification_base)
        oligo.build()



class edit_json_dialog():
    def __init__(self, init_json, input_label):
        if init_json == '':
            self.json = json.dumps({'key 1': 'value 1'})
        else:
            self.json = init_json
        self.input_label = input_label
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                ui.input(label='json type', value=self.input_label)
                self.editor = ui.json_editor({'content': {'json': self.json}}).style(
                    'width: 600px; height: 800px; font-size: 20px')
                with ui.row():
                    ui.button('Save', on_click=self.get_data)
                    ui.button('Close', on_click=self.do_close)


    def on_save(self, data):
        ui.notify(data)

    async def get_data(self):
        data = await self.editor.run_editor_method('get')
        self.on_save(data)


    def do_close(self):
        self.dialog.close()


if __name__ in {"__main__", "__mp_main__"}:
    run_mods()