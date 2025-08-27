from nicegui import events, ui, app
from OligoMap_utils import api_db_interface
import requests
import pandas as pd
import json
from pubchempy import get_compounds, Compound, get_substances, get_properties, get_synonyms

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



class chemicals_page_model(api_db_interface):
    def __init__(self, api_IP, db_port):
        super().__init__(api_IP, db_port)

        if 'pincode' in list(app.storage.user.keys()):
            self.pincode = app.storage.user.get('pincode')
        else:
            self.pincode = ''

        self.pubchem = pubchem_interface()