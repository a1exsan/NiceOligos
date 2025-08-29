from nicegui import events, ui, app
from OligoMap_utils import api_db_interface
import requests
import pandas as pd
import json
from pubchempy import get_compounds, Compound, get_substances, get_properties, get_synonyms
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, rdMolDescriptors
import io
import base64
from pprint import pformat
from rdkit.Chem.Draw import rdMolDraw2D

class pubchem_interface():
    def __init__(self):
        self.draw()

    def draw(self):
        with ui.card():
            with ui.row():
                self.search_text = ui.input('Chemical name',
                        placeholder='Enter chemical name').classes('w-[400px]').style('font-size: 20px;')
                self.search_button = ui.button('Search', on_click=self.on_search_click).classes('w-[200px]')

            self.search_result = ui.textarea().classes('w-full h-full large-textarea')

    def search(self, text, params={'namespace': 'smiles', 'searchtype': 'smiles'}):
        #return get_compounds(text, namespace=params['namespace'], searchtype=params['searchtype'])
        return get_compounds(text, "name")
        #return get_properties(text, "name")

    def on_search_click(self):
        text = self.search_text.value
        #comps = self.search(text)

        syns = get_synonyms(text, "name")

        d = {'cid': syns[0]['CID']}
        d['Synonym'] = syns[0]['Synonym']
        vioxx = Compound.from_cid(d['cid'])
        d['molecular_formula'] = vioxx.molecular_formula
        d['molecular_weight'] = vioxx.molecular_weight
        d['inchi'] = vioxx.inchi
        d['xlogp'] = vioxx.xlogp

        out = f'{text}: {d}'

        self.search_result.value = out

        #ui.notify(out)

class phys_chem_props_interface():
    def __init__(self):
        self.draw()

    def draw(self):
        with ui.card():

            with ui.column():
                self.compound_name = ui.input(label='Название материала:').classes('w-[400px]').style('font-size: 20px;')
                self.producer = ui.input(label='Производитель:').classes('w-[400px]').style('font-size: 20px;')
                self.supplyer = ui.input(label='Поставщик:').classes('w-[400px]').style('font-size: 20px;')
                with ui.row():
                    self.price = ui.input(label='Цена за единицу:').classes('w-[400px]').style('font-size: 20px;')
                    self.units = ui.input(label='Единица измерения:').classes('w-[400px]').style('font-size: 20px;')
                self.description = ui.textarea(label='Общее описание:').classes('w-[1000px]').style('font-size: 18px;')

            ui.link('Открыть редактор структур Ketcher',
                    'https://lifescience.opensource.epam.com/KetcherDemoSA/index.html', new_tab=True
                    ).style('font-size: 20px;')
            ui.label('Введите тектовое описание структуры (InChi)').style('font-size: 20px;')
            self.structure = ui.textarea(on_change=self.draw_structure).classes('w-[1000px]')

            with ui.row():
                self.props = {}
                self.mol_svg = ui.html('')
                #self.mol_image = ui.image('').style('width: 600px; height: 400px; margin-top: 20px;')
                with ui.column():
                    self.mol_props = ui.textarea(label='Свойства').classes('w-[400px]').style('font-size: 18px;')
                    self.mol_description = ui.textarea(label='Введите описание от Lumiprobe'
                        #on_change=self.on_change_mol_description
                    ).classes('w-[400px]').style('font-size: 18px;')
                    ui.button('set', on_click=self.parse_mol_description)

            self.structure_adduct = ui.textarea(label='Введити InChi аддукта (это данные для массспектрометрии)',
                                                on_change=self.draw_adduct_structure).classes('w-[1000px]')
            with ui.row():
                self.adduct_props = {}
                self.mol_adduct_svg = ui.html('')
                #self.mol_image = ui.image('').style('width: 600px; height: 400px; margin-top: 20px;')
                with ui.column():
                    self.mol_adduct_props = ui.textarea(label='Свойства аддукта').classes('w-[400px]').style('font-size: 18px;')

            with ui.row():
                self.save_data = ui.button('Сохранить', color='green',
                                       on_click=self.save_product_data_to_base).classes('w-[200px]')
                self.save_data = ui.button('Отмена', color='orange',
                                       on_click=self.on_cencel).classes('w-[200px]')


    def on_save(self, data):
        pass

    def save_product_data_to_base(self):
        desc = {}
        desc['name'] = self.compound_name.value
        desc['producer'] = self.producer.value
        desc['supplyer'] = self.supplyer.value
        desc['price'] = self.price.value
        desc['price_units'] = self.units.value
        desc['price_date'] = datetime.now().date().strftime('%d.%m.%Y')
        desc['price_date_format'] = '%d.%m.%Y'
        desc['mol_inchi'] = self.structure.value
        desc['mol_props'] = json.dumps(self.props)
        desc['mol_lumiprobe_data'] = json.dumps(self.mol_lumiprobe_dict)
        desc['adduct_inchi'] = self.structure_adduct.value
        desc['adduct_props'] = json.dumps(self.adduct_props)

        self.on_save(desc)
        #print(json.dumps(desc))

    def on_cencel(self):
        pass

    def parse_mol_description(self):
        s = self.mol_description.value
        if s != '':
            d = {}
            for row in s.split('\n'):
                d[row[:row.find(':')]] = row[row.find(':') + 1:].replace('\t', '')
            self.mol_description.value = str(d)
            self.mol_lumiprobe_dict = d.copy()
            #print(self.mol_lumiprobe_dict)

    def culc_structure(self, inchi):
        try:
            #img_data = self.draw_mol_from_inchi(self.structure.value)
            #self.mol_image.set_source(img_data)
            props = self.get_physchem_properties(inchi)
            svg = self.draw_mol_svg(inchi)
            return props, svg
        except Exception as e:
            ui.notify(f'Ошибка: {e}', color='red')

    def draw_structure(self):
        self.props, self.mol_svg.content = self.culc_structure(self.structure.value)
        self.mol_props.value = str(self.props)

    def draw_adduct_structure(self):
        self.adduct_props, self.mol_adduct_svg.content = self.culc_structure(self.structure_adduct.value)
        self.mol_adduct_props.value = str(self.adduct_props)

    def draw_mol_from_inchi(self, inchi):
        mol = Chem.inchi.MolFromInchi(inchi)
        if mol is None:
            raise ValueError("Некорректный InChI")
        img = Draw.MolToImage(mol, size=(600, 400))
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    def draw_mol_svg(self, inchi):
        mol = Chem.inchi.MolFromInchi(inchi)
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 400)
        rdMolDraw2D.SetDarkMode(drawer)
        #opts = drawer.drawOptions()
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()


    def get_physchem_properties(self, inchi):
        mol = Chem.inchi.MolFromInchi(inchi)
        if mol is None:
            raise ValueError('Некорректный InChI')
        return {
            'Mol weight, Da ': round(Descriptors.MolWt(mol), 2),
            'Brutto ': rdMolDescriptors.CalcMolFormula(mol),
            #'ЛогP (октанол/вода)': Descriptors.MolLogP(mol),
            #'Доноры водородных связей': rdMolDescriptors.CalcNumHBD(mol),
            #'Акцепторы водородных связей': rdMolDescriptors.CalcNumHBA(mol),
            #'TPSA (полярная поверхность)': rdMolDescriptors.CalcTPSA(mol),
            #'Число ротируемых связей': Descriptors.NumRotatableBonds(mol),
            #'Число колец': Descriptors.RingCount(mol),
        }



class chemicals_page_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        else:
            self.pincode = ''

        self.pubchem = pubchem_interface()
        self.pc_props = phys_chem_props_interface()