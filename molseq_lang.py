
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
import pandas as pd
import random
import threading
import time
from synthesis_method import synth_base


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


    def draw_context(self, context, width=600, height=400):
        context.content = ''
        context.content += (f'<rect x={0} y={0} '
                            f'width={self.context_width} height={self.context_height} '
                            f' rx=10 ry=10 fill="{None}" fill-opacity="{0.6}"'
                            f'stroke="{self.border_color}" stroke-width="4"/>')
        if self.rnx_sel_id > 0:
            self.reaction_base[self.rnx_sel_id].draw_rnx_svg(context, width=width, height=height)
        if self.mod_sel_id > 0:
            self.modification_base[self.mod_sel_id].draw_mod_svg(context, width=width, height=height)
            if self.reactant.smiles != '':
                self.reactant.draw_mod_svg_reactant(context, width=width, height=height)

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

    def del_mod_from_chain(self, mod_symbol):
        if mod_symbol in self.chain:
            self.chain.remove(mod_symbol)
        return ''.join(self.chain)

class single_nucleic_acid_chain_assembler(single_nucleic_acid_chain):
    def __init__(self, seq, rnx_base, mod_base):
        super().__init__(seq)
        self.do_when_build_finished = None
        self.build_progress = None
        self.rnx_base = rnx_base
        self.DMT_on = False
        base = {}
        rowdata = []
        for val in mod_base.values():
            base[val.symbol] = val
            rowdata.append(val.to_dict())
        self.mod_base = base
        self.structure = ''
        self.scheme = {}
        self.desupport_id = 0
        self.border_color = 'gray'
        self.mod_base_df = pd.DataFrame(rowdata)
        df = self.mod_base_df[self.mod_base_df['symbol'].str.contains('_class')]
        self.mod_classes = df.to_dict('records')

    def check_chain(self):
        errors = {}
        class_chain = [json.loads(self.mod_base[mod].data_json)['class'] for mod in self.chain[::-1]]
        if class_chain[0] == 'CPG' and class_chain[1] == 'amidite':
            index = 1
            amidite_end = False
            while True:
                if amidite_end:
                    if class_chain[index] == 'amidite':
                        errors['5end'] = 'missing amidite end'
                if class_chain[index] != 'amidite':
                    amidite_end = True
                index += 1
                if index == len(class_chain):
                    break
        else:
            if class_chain[0] == 'NHS':
                errors['3end'] = 'NHS'
            elif class_chain[0] == 'azide':
                errors['3end'] = 'azide'
            elif class_chain[0] == 'amidite':
                errors['3end'] = 'no CPG'
        return json.dumps(errors)

    def draw_structure(self, context, width, height):
        context.content = ''
        context.content += (f'<rect x={0} y={0} '
                            f'width={width} height={height} '
                            f' rx=10 ry=10 fill="{None}" fill-opacity="{0.6}"'
                            f'stroke="{self.border_color}" stroke-width="4"/>')
        mol = moleculeInfo(self.structure)
        svg = mol.draw_svg(width-20, height-20)
        context.content += f'<g transform="translate(10, 10)">{svg}</g>'

    def get_branch_count(self, structure):
        molecule = Chem.MolFromSmiles(structure)
        substructure = Chem.MolFromSmiles('C1C=C(C(O)(C2C=CC=CC=2)C2C=CC=CC=2)C=CC=1')
        matches = molecule.GetSubstructMatches(substructure)
        if len(matches) == 0:
            return 1
        else:
            return len(matches)

    def get_structure_count(self, structure, substructure):
        molecule = Chem.MolFromSmiles(structure)
        substructure_ = Chem.MolFromSmiles(substructure)
        matches = molecule.GetSubstructMatches(substructure_)
        if len(matches) == 0:
            return 1
        else:
            return len(matches)

    def do_auto_cycle_detrit(self, mod_data):
        if 'DETRIT' in mod_data:
            for id in mod_data['DETRIT']:
                rnx_smarts = self.rnx_base[id].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_auto_cycle_couple(self, mod_data, adduct):
        if 'COUPLE' in mod_data:
            for id in mod_data['COUPLE']:
                rnx_smarts = self.rnx_base[id].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                react2 = Chem.MolFromSmiles(self.structure)
                react1 = Chem.MolFromSmiles(adduct)
                products = rxn.RunReactants([react1, react2])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_auto_cycle_oxid(self, mod_data):
        if 'OXID' in mod_data:
            for id in mod_data['OXID']:
                rnx_smarts = self.rnx_base[id].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                ox = Chem.MolFromSmiles('O')
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react, ox])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_auto_cycle_debl(self, mod_data):
        if 'DEBL' in mod_data:
            for id in mod_data['DEBL']:
                rnx_smarts = self.rnx_base[id].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_auto_reactions(self, modification_):
        if self.structure == '':
            self.structure = modification_.smiles
        smiles_list = modification_.smiles.split('.')
        adduct = smiles_list[0]
        mod_data = json.loads(modification_.data_json)
        branch_count = self.get_branch_count(self.structure)
        if 'DESUPPORT' in mod_data:
            self.desupport_id = mod_data['DESUPPORT'][0]
        #print(modification_.to_dict())
        #print(branch_count)
        for i in range(branch_count):
            self.do_auto_cycle_detrit(mod_data=mod_data)
            #print(self.structure)
        for i in range(branch_count):
            self.do_auto_cycle_couple(mod_data=mod_data, adduct=adduct)
            #print(self.structure)
        for i in range(branch_count):
            self.do_auto_cycle_oxid(mod_data=mod_data)
        for i in range(branch_count):
            self.do_auto_cycle_debl(mod_data=mod_data)


    def do_desupport_structure(self):
        if self.structure != '':
            if self.desupport_id > 0:
                rnx_smarts = self.rnx_base[self.desupport_id].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_final_detrit_structure(self, DMT_on=True):
        if self.structure != '':
            if not DMT_on:
                rnx_smarts = self.rnx_base[1].smarts
                rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                react = Chem.MolFromSmiles(self.structure)
                products = rxn.RunReactants([react])
                if products != ():
                    self.structure = Chem.MolToSmiles(products[0][0])

    def do_click_reaction_on_structure(self, modification_):
        if self.structure != '':
            mod_data = json.loads(modification_.data_json)
            if 'click' in mod_data:
                repeats = 1
                if 'class' in mod_data:
                    if mod_data['class'] == 'NHS':
                        repeats = self.get_structure_count(self.structure, 'CCN')
                    elif mod_data['class'] == 'azide':
                        repeats = self.get_structure_count(self.structure, 'C#CC')
                for i in range(repeats):
                    for id in mod_data['click']:
                        rnx_smarts = self.rnx_base[id].smarts
                        rxn = rdChemReactions.ReactionFromSmarts(rnx_smarts)
                        react1 = Chem.MolFromSmiles(modification_.smiles)
                        react2 = Chem.MolFromSmiles(self.structure)
                        products = rxn.RunReactants([react1, react2])
                        if products != ():
                            self.structure = Chem.MolToSmiles(products[0][0])
                        else:
                            products = rxn.RunReactants([react2, react1])
                            if products != ():
                                self.structure = Chem.MolToSmiles(products[0][0])

    def build(self):
        reverse_chain = self.chain[::-1]
        self.structure = ''
        self.desupport_id = 0

        if self.build_progress is not None:
            self.build_progress.value = 0
            time.sleep(0.5)

        for i, token in enumerate(reverse_chain):
            if self.build_progress is not None:
                self.build_progress.value = round(i / len(reverse_chain), 2)
                time.sleep(0.1)
            if token in self.mod_base:
                #print(token, self.get_structure_class(self.mod_base[token].smiles))
                self.do_auto_reactions(self.mod_base[token])
                self.do_click_reaction_on_structure(self.mod_base[token])
            else:
                print(f'{token} not in base')
        self.do_desupport_structure()
        self.do_final_detrit_structure(DMT_on=self.DMT_on)

        if self.build_progress is not None:
            self.build_progress.value = 1
            time.sleep(0.1)

        if self.do_when_build_finished is not None:
            self.do_when_build_finished({'structure': self.structure})

    def run_build(self):
        background_thread = threading.Thread(target=self.build)
        background_thread.daemon = True
        background_thread.start()


    def get_structure_class(self, smiles):
        result = ''
        if smiles != '':
            molecule = Chem.MolFromSmiles(smiles)
            for row in self.mod_classes:
                substructure = Chem.MolFromSmiles(row['smiles'])
                if molecule.HasSubstructMatch(substructure):
                    matches = molecule.GetSubstructMatches(substructure)
                    data = json.loads(row['data_json'])
                    if 'class' in data:
                        result = data['class'], len(matches)
                        break
        return result

    def compile_structure(self):
        reverse_chain = self.chain[::-1]
        self.scheme = {}
        self.scheme['unknown_symbol'] = []
        self.scheme['unknown_class'] = []
        self.scheme['error_3end'] = []
        self.scheme['rowdata'] = []
        self.scheme['auto_sequence'] = None
        self.scheme['modif_overlimit'] = None
        for i, token in enumerate(reverse_chain):
            if token not in self.mod_base:
                self.scheme['unknown_symbol'].append(token)
            else:
                j_data = json.loads(self.mod_base[token].data_json)
                if 'class' in j_data:
                    t_class = j_data['class']
                    d = {}
                    d['symbol'] = token
                    d['class'] = t_class
                    d['index'] = i
                    self.scheme['rowdata'].append(d)
                    if i == 0 and t_class != 'CPG':
                        self.scheme['error_3end'].append(token)
                else:
                    self.scheme['unknown_class'].append(token)
        if len(self.scheme['unknown_symbol']) == 0 and len(self.scheme['unknown_class']) == 0 and len(self.scheme['error_3end'])==0:
            df = pd.DataFrame(self.scheme['rowdata'])
            if df.shape[0] > 0:
                classes_set = list(df['class'].unique())
                if classes_set == ['amidite', 'CPG'] or classes_set == ['CPG', 'amidite']:
                    self.scheme['auto_sequence'] = self.chain
                elif len(classes_set) <= 5:
                    df.sort_values(by='index', ascending=False, inplace=True)
                    self.scheme['auto_sequence'] = self.compile_auto_sequence(df)
                else:
                    self.scheme['modif_overlimit'] = len(classes_set)

    def compile_auto_sequence(self, dataframe):
        cls_l, symb_l = list(dataframe['class']), list(dataframe['symbol'])
        seq = []
        for i in range(dataframe.shape[0]):
            if cls_l[i] == 'azide' and i == 0:
                if cls_l[i+1] == 'amidite' or cls_l[i+1] == 'CPG' or cls_l[i+1] == 'NHS' or cls_l[i+1] == 'alkine':
                    seq.append('[Alk]')
                elif cls_l[i+1] == 'azide':
                    seq.append('[Alk_dmt]')
                else:
                    seq.append('[Alk_dmt]')
            elif cls_l[i] == 'azide' and i > 0:
                if cls_l[i - 1] == 'azide':
                    seq.append('[Alk_dmt]')
                else:
                    seq.append('[Alk_dT]')
            #elif cls_l[i] == 'azide' and i == dataframe.shape[0] - 1:
            #    seq.append('[NH2_cpg500]')
            if cls_l[i] == 'NHS' and i == 0:
                if cls_l[i + 1] == 'amidite' or cls_l[i + 1] == 'CPG' or cls_l[i + 1] == 'azide' or cls_l[i + 1] == 'alkine':
                    seq.append('[NH2_C6]')
                elif cls_l[i + 1] == 'NHS':
                    seq.append('[Alk_dmt]')
                else:
                    seq.append('[NH2_C6]')
            elif cls_l[i] == 'NHS' and i > 0:
                if cls_l[i - 1] == 'NHS':
                    seq.append('[Alk_dmt]')
                else:
                    seq.append('[Alk_dT]')
            elif cls_l[i] == 'NHS' and i == dataframe.shape[0] - 1:
                seq.append('[NH2_cpg500]')
            if cls_l[i] in ['amidite', 'CPG']:
                seq.append(symb_l[i])
        return seq





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
            {"field": "symbol", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "smiles", 'editable': True,  'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "unicode", 'editable': True,  'filter': 'agTextColumnFilter', 'floatingFilter': True},
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
            with ui.column():
                ui.button('Parce seq', on_click=self.on_parse_seq)
                self.check_dmt = ui.checkbox('DMT on', value=True)
                ui.button('Compile', on_click=self.on_compile_seq)
                ui.button('Edit sheme', on_click=self.on_edit_synth_scheme)
            self.seq_area = ui.textarea(label='Sequence').style('width: 400px; font-size: 18px')
            with ui.column():
                self.tokens_area = ui.textarea(label='Product smiles').style('width: 800px')
                self.progress = ui.linear_progress().style('width: 800px')
            self.props_area = ui.textarea(label='Mol props').style('width: 400px; font-size: 20px;')

        self.react_width = 2400
        self.react_height = 1600
        bkg = image_background(self.react_width, self.react_height, color=self.bkg_color)
        self.react_image = ui.interactive_image(f"data:image/png;base64,{bkg.background_base64}",
                                              size=(self.react_width, self.react_height), on_mouse=self.mouse_handler,
                                              events=['mousedown', 'mousemove', 'mouseup', 'click']
                                              )
        self.react_context = self.react_image.add_layer()

        self.init_reaction_rowdata()

        self.method_base = synth_base()
        with ui.row():
            colDefs = [
                {"field": "id", 'editable': False},
                {"field": "synth_name", 'editable': True},
                {"field": "scale", 'editable': True},
                {"field": "data_json", 'editable': True},
                {"field": "template", 'editable': True}
            ]
            with ui.column():
                ui.button('Add method', on_click=self.on_add_synth_method)
                ui.button('load method', on_click=self.on_load_synth_method)
                ui.button('save method', on_click=self.on_save_synth_method)
            self.method_grid = ui.aggrid(
                {
                    'columnDefs': colDefs,
                    'rowData': [],
                    'rowSelection': 'multiple',
                    "pagination": True,
                    "enterNavigatesVertically": True,
                    "enterNavigatesVerticallyAfterEdit": True,
                    "singleClickEdit": True,
                    # "enableRangeSelection": True,
                },
                theme='alpine-dark').style('height: 600px; width: 1200px')

            self.method_grid.auto_size_columns = True
            self.method_grid.on("cellValueChanged", self.update_grid_cell_data_method)
            self.method_grid.on('rowSelected', self.on_select_method_row)
            self.method_grid.options['rowData'] = self.method_base.get_all_rowdata()
            self.method_grid.update()

            self.method_json = {}
            self.method_jedit = ui.json_editor({'content': {'json': self.method_json}}).style(
                'width: 600px; height: 600px; font-size: 20px')



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

    def on_add_synth_method(self):
        data = {'synth_name': 'default',
             'scale': '3 mg',
             'data_json': '',
             'template': ''}
        self.method_base.insert_method(data)
        self.method_grid.options['rowData'] = self.method_base.get_all_rowdata()
        self.method_grid.update()

    def on_load_synth_method(self):
        if self.method_base.selected_id > 0:
            _json = self.method_grid.options['rowData'][self.method_base.selected_id - 1]['data_json']
            json_obj = json.loads(_json)
            if isinstance(json_obj, dict) and 'text' in json_obj:
                json_obj = json.loads(json_obj['text'])
            self.method_json = json_obj
            #self.method_jedit.set_value(json_obj)
            print(self.method_json)

            self.method_jedit.content['json'] = json_obj

            self.method_jedit.update()


    async def on_save_synth_method(self):
        if self.method_base.selected_id > 0:
            d = self.method_grid.options['rowData'][self.method_base.selected_id - 1]
            jd = await self.method_jedit.run_editor_method('get')
            data = {'synth_name': d['synth_name'],
                'scale': d['scale'],
                'data_json': json.dumps(jd),
                'template': d['template']}
            self.method_base.update_method(self.method_base.selected_id, data)


    def update_grid_cell_data_method(self, e):
        data = {'synth_name': e.args['data']['synth_name'],
             'scale': e.args['data']['scale'],
             'data_json': e.args['data']['data_json'],
             'template': e.args['data']['template']}
        self.method_base.update_method(e.args['data']['id'], data)

    async def on_select_method_row(self):
        selrows = await self.method_grid.get_selected_rows()
        self.method_base.selected_id = selrows[0]['id']
        ui.notify(f'Выбран метод: {selrows[0]["synth_name"]} масштаб {selrows[0]["scale"]}')


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
            oligo = single_nucleic_acid_chain_assembler('ACGT',
                                                        self.obj_base.reaction_base,
                                                        self.obj_base.modification_base)
            data_json = json.loads(data['json'])
            data_json['class'], data_json['class_count'] = oligo.get_structure_class(mod_dict['smiles'])
            d = {'symbol': mod_dict['symbol'],
                 'unicode': mod_dict['unicode'],
                 'smiles': mod_dict['smiles'],
                 'data_json': json.dumps(data_json)}
            mod = modification.from_dict(d)
            self.obj_base.update_modification(self.obj_base.mod_sel_id, mod)
            self.init_reaction_rowdata()

    def on_build_finished(self, data):
        self.tokens_area.value = data['structure']
        oligo = single_nucleic_acid_chain_assembler('ACGT',
                                                    self.obj_base.reaction_base,
                                                    self.obj_base.modification_base)
        oligo.structure = data['structure']
        oligo.draw_structure(self.react_context, self.react_width, self.react_height)
        mol = moleculeInfo(data['structure'])
        self.props_area.value = json.dumps(mol.get_props())

    def on_parse_seq(self):
        sequence = self.seq_area.value
        oligo = single_nucleic_acid_chain(sequence)
        self.tokens_area.value = json.dumps(oligo.chain)

        oligo = single_nucleic_acid_chain_assembler(sequence,
                                                      self.obj_base.reaction_base,
                                                      self.obj_base.modification_base)
        oligo.DMT_on = self.check_dmt.value
        oligo.build_progress = self.progress
        oligo.do_when_build_finished = self.on_build_finished
        oligo.run_build()


    def on_compile_seq(self):
        sequence = self.seq_area.value
        oligo = single_nucleic_acid_chain(sequence)
        self.tokens_area.value = json.dumps(oligo.chain)

        oligo = single_nucleic_acid_chain_assembler(sequence,
                                                      self.obj_base.reaction_base,
                                                      self.obj_base.modification_base)
        oligo.compile_structure()
        self.tokens_area.value = json.dumps(oligo.scheme)

    def on_edit_synth_scheme(self):
        rowdata = []
        dl = 'A C G T N R Y'.split(' ')
        for i in range(30):
            d = {}
            d['#'] = i + 1
            d['Sequence'] = ''.join([dl[random.randint(0,6)] for i in range(random.randint(20,50))])
            if random.randint(1,5) == 1:
                d['Sequence'] = d['Sequence'] + '[BHQ1]'
            elif random.randint(1,5) == 2:
                d['Sequence'] = d['Sequence'] + '[BHQ3]'
            d['Position'] = 'A1'
            rowdata.append(d)
        sheme = synth_scheme_dialog(rowdata, self.obj_base)
        sheme.dialog.open()

class reagent_tab():
    def __init__(self, rowdata, mod_base):
        self.rowdata = rowdata
        self.mod_base = {}
        for val in mod_base.values():
            self.mod_base[val.symbol] = val

    def get_wobblw_count(self, w):
        A, C, G, T = 'A', 'C', 'G', 'T'
        if w == 'R':
            return [A, G]
        if w == 'Y':
            return [C, T]
        if w == 'K':
            return [T, G]
        if w == 'M':
            return [A, C]
        if w == 'S':
            return [C, G]
        if w == 'W':
            return [A, T]
        if w == 'B':
            return [C, G, T]
        if w == 'D':
            return [A, G, T]
        if w == 'H':
            return [A, C, T]
        if w == 'V':
            return [A, C, G]
        if w == 'N':
            return [A, C, G, T]

    def get_reagents(self):
        out = []
        d = {}
        for row in self.rowdata:
            oligo = single_nucleic_acid_chain(row['Chain'])
            for cell in oligo.chain:
                if cell not in d:
                    d[cell] = 0
                else:
                    d[cell] += 1
        for key, value in zip(d.keys(), d.values()):
            row = {}
            row['symbol'] = key
            row['count'] = value
            if key in self.mod_base:
                row['unicode'] = self.mod_base[key].unicode
            else:
                row['unicode'] = ''
            out.append(row)

        return out


class synth_scheme_dialog():
    def __init__(self, rowdata, obj_base):
        self.rowdata = rowdata
        self.obj_base = obj_base
        self.context_width = 800
        self.context_height = 400
        self.init_rowdata = []
        self.get_init_rowdata()
        self.init_dialog()

    def get_init_rowdata(self):
        self.init_rowdata = []
        for row in self.rowdata:
            d = {}
            d['#'] = row['#']
            d['Position'] = row['Position']
            d['Sequence'] = row['Sequence']
            d['DMT on'] = True
            d['ASM sequence'] = ''
            d['Chain'] = d['Sequence']
            oligo = single_nucleic_acid_chain_assembler(d['Chain'],
                                                        self.obj_base.reaction_base,
                                                        self.obj_base.modification_base)
            d['errors'] = json.loads(oligo.check_chain())
            if d['errors'] != {}:
                if '3end' in d['errors']:
                    if d['errors']['3end'] == 'no CPG':
                        d['Chain'] = d['Chain'] + self.get_cpg_mod(d['Chain'])
                        oligo.chain = oligo.parse(d['Chain'])
            d['errors'] = oligo.check_chain()
            self.init_rowdata.append(d)

    def get_cpg_mod(self, chain):
        if len(chain) <= 30:
            return '[CPG500]'
        elif len(chain) > 30 and len(chain) <= 65:
            return '[CPG1000]'
        elif len(chain) > 65 and len(chain) <= 110:
            return '[CPG2000]'
        elif len(chain) > 110 and len(chain) <= 500:
            return '[CPG3000]'
        else:
            return '[CPG3000]'

    def init_dialog(self):
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                with ui.row():
                    ui.button("add 5' mod", color="green", on_click=self.on_add_5mod)
                    ui.button("add 3' mod", color="green", on_click=self.on_add_3mod)
                    ui.button("del mod", color="red", on_click=self.on_del_mod)
                    ui.button("on dmt", color="green", on_click=self.on_dmt_on)
                    ui.button("off dmt", color="red", on_click=self.on_dmt_off)
                    ui.button("get reagents", on_click=self.on_get_reagent_tab)
                    colDefs = [
                        {"field": "#", 'editable': False},
                        {"field": "Position", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "Sequence", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "DMT on", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "Chain", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "ASM sequence", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                        {"field": "errors", 'editable': False},
                    ]
                    colDefs_reagent = [
                        {"field": "symbol", 'editable': True},
                        {"field": "count", 'editable': False},
                        {"field": "asm reagent", 'editable': True},
                        {"field": "Conc, mg/ml", 'editable': True},
                        {"field": "Reagent, mg", 'editable': True},
                        {"field": "Reagent, ml", 'editable': True},
                        {"field": "Consumption", 'editable': True},
                        {"field": "unicode", 'editable': False},
                    ]
                with ui.row():
                        self.scheme_grid = ui.aggrid(
                        {
                            'columnDefs': colDefs,
                            'rowData': self.init_rowdata,
                            'rowSelection': 'multiple',
                            "pagination": True,
                            "enterNavigatesVertically": True,
                            "enterNavigatesVerticallyAfterEdit": True,
                            "singleClickEdit": True,
                            # "enableRangeSelection": True,
                        },
                        theme='alpine-dark').style('height: 900px; width: 1600px')  # alpine  material  quartz  balham
                        self.scheme_grid.auto_size_columns = True
                        self.scheme_grid.on("cellValueChanged", self.update_grid_cell_data_scheme)

                        self.reagent_grid = ui.aggrid(
                            {
                                'columnDefs': colDefs_reagent,
                                'rowData': [],
                                'rowSelection': 'multiple',
                                "pagination": True,
                                "enterNavigatesVertically": True,
                                "enterNavigatesVerticallyAfterEdit": True,
                                "singleClickEdit": True,
                                # "enableRangeSelection": True,
                            },
                            theme='alpine-dark').style(
                            'height: 900px; width: 900px')  # alpine  material  quartz  balham
                        self.reagent_grid.auto_size_columns = True
                        self.reagent_grid.on("cellValueChanged", self.update_reagent_tab)

                with ui.column():
                    with ui.row():
                        mod_rowdata = self.obj_base.get_modification_rowdata()
                        colDefs_modif = [
                            {"field": "id", 'editable': False},
                            {"field": "symbol", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                            {"field": "smiles", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                            {"field": "unicode", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
                            {"field": "data_json", 'editable': True}
                        ]
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
                            theme='alpine-dark').style('height: 400px; width: 800px')
                        self.mod_grid.auto_size_columns = True
                        self.mod_grid.on('rowSelected', self.on_select_mod_row)

                        self.img_width = self.context_width
                        self.img_height = self.context_height
                        self.bkg_color = '#010101'
                        self.border_color = 'gray'
                        bkg = image_background(self.img_width, self.img_height, color=self.bkg_color)
                        self.image = ui.interactive_image(f"data:image/png;base64,{bkg.background_base64}",
                                                          size=(self.img_width, self.img_height),
                                                          on_mouse=self.mouse_handler,
                                                          events=['mousedown', 'mousemove', 'mouseup', 'click']
                                                          )
                        self.rnx_context = self.image.add_layer()

                with ui.row():
                    ui.button('Save', on_click=self.get_data)
                    ui.button('Close', on_click=self.do_close)

    def on_save(self, data):
        ui.notify(data)

    def get_data(self):
        self.on_save('')

    def do_close(self):
        self.dialog.close()

    async def on_add_5mod(self):
        mod_row = await self.mod_grid.get_selected_rows()
        scheme_row = await self.scheme_grid.get_selected_rows()
        df = pd.DataFrame(scheme_row)
        sel_list = list(df['#'])
        rowdata = self.scheme_grid.options['rowData']
        for i, row in enumerate(rowdata):
            if row['#'] in sel_list:
                d = row.copy()
                d['Chain'] = mod_row[0]['symbol'] + row['Chain']
                self.scheme_grid.options['rowData'][i] = d
        self.scheme_grid.update()

    async def on_del_mod(self):
        mod_row = await self.mod_grid.get_selected_rows()
        mod = mod_row[0]['symbol']
        scheme_row = await self.scheme_grid.get_selected_rows()
        df = pd.DataFrame(scheme_row)
        sel_list = list(df['#'])
        rowdata = self.scheme_grid.options['rowData']
        for i, row in enumerate(rowdata):
            if row['#'] in sel_list:
                d = row.copy()
                oligo = single_nucleic_acid_chain(row['Chain'])
                d['Chain'] = oligo.del_mod_from_chain(mod)
                self.scheme_grid.options['rowData'][i] = d
        self.scheme_grid.update()


    async def on_add_3mod(self):
        mod_row = await self.mod_grid.get_selected_rows()
        scheme_row = await self.scheme_grid.get_selected_rows()
        df = pd.DataFrame(scheme_row)
        sel_list = list(df['#'])
        rowdata = self.scheme_grid.options['rowData']
        for i, row in enumerate(rowdata):
            if row['#'] in sel_list:
                d = row.copy()
                d['Chain'] = row['Chain'] + mod_row[0]['symbol']
                self.scheme_grid.options['rowData'][i] = d
        self.scheme_grid.update()


    async def on_dmt_on(self):
        scheme_row = await self.scheme_grid.get_selected_rows()
        df = pd.DataFrame(scheme_row)
        sel_list = list(df['#'])
        rowdata = self.scheme_grid.options['rowData']
        for i, row in enumerate(rowdata):
            if row['#'] in sel_list:
                d = row.copy()
                d['DMT on'] = True
                self.scheme_grid.options['rowData'][i] = d
        self.scheme_grid.update()


    async def on_dmt_off(self):
        scheme_row = await self.scheme_grid.get_selected_rows()
        df = pd.DataFrame(scheme_row)
        sel_list = list(df['#'])
        rowdata = self.scheme_grid.options['rowData']
        for i, row in enumerate(rowdata):
            if row['#'] in sel_list:
                d = row.copy()
                d['DMT on'] = False
                self.scheme_grid.options['rowData'][i] = d
        self.scheme_grid.update()


    def on_get_reagent_tab(self):
        tab = reagent_tab(self.scheme_grid.options['rowData'], mod_base=self.obj_base.modification_base)
        self.reagent_grid.options['rowData'] = tab.get_reagents()
        self.reagent_grid.update()


    def mouse_handler(self, e):
        pass

    def draw_context(self, context, id):
        context.content = ''
        context.content += (f'<rect x={0} y={0} '
                            f'width={self.img_width} height={self.img_height} '
                            f' rx=10 ry=10 fill="{None}" fill-opacity="{0.6}"'
                            f'stroke="{self.border_color}" stroke-width="4"/>')
        mol = moleculeInfo(self.obj_base.modification_base[id].smiles)
        svg = mol.draw_svg(self.img_width-20, self.img_height-20)
        context.content += f'<g transform="translate(10, 10)">{svg}</g>'

    async def on_select_mod_row(self):
        selrows = await self.mod_grid.get_selected_rows()
        self.draw_context(self.rnx_context, selrows[0]['id'])

    def update_grid_cell_data_scheme(self, e):
        self.scheme_grid.options['rowData'][e.args['rowIndex']] = e.args['data']
        self.scheme_grid.update()

    def update_reagent_tab(self, e):
        self.reagent_grid.options['rowData'][e.args['rowIndex']] = e.args['data']
        self.reagent_grid.update()



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