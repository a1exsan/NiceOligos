import re
from chemicals_page import moleculeInfo
from OligoMap_utils import api_db_interface
import json
import requests

class modification():
    def __init__(self, symbol, smiles=None, unicode=None):
        self.symbol = symbol
        self.smiles = smiles
        self.unicode = unicode

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



class modification_base(api_db_interface):
    def __init__(self, db_IP, db_port):
        super().__init__(db_IP, db_port)
        self._base = {}

    def insert_modification(self, data_object):
        url = f'{self.api_db_url}/insert_data/{self.mod_db_name}/mod'
        r = requests.post(url,
                          json=json.dumps([
                              data_object['symbol'],
                              data_object['unicode'],
                              data_object['smiles'],
                              data_object['data_json']
                          ]), headers=self.headers())

    def insert_reaction(self, data_object):
        url = f'{self.api_db_url}/insert_data/{self.mod_db_name}/rnx'
        r = requests.post(url,
                          json=json.dumps([
                              data_object['name'],
                              data_object['smarts'],
                              data_object['data_json']
                          ]), headers=self.headers())


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
    mod_base.insert_reaction(rnx_obj)

    rnx_obj = {}
    rnx_obj['name'] = 'COUPLING main'
    rnx_obj[
        'smarts'] = '[#8:1]-[P@;H0;D3;+0:2](-[#8:3])-[N;H0;D3;+0:4](-[C:5])-[C:6].[C:7]-[OH;D1;+0:8]>>[#8:1]-[P;H0;D3;+0:2](-[#8:3])-[O;H0;D2;+0:8]-[C:7].[C:5]-[NH;D2;+0:4]-[C:6]'
    rnx_obj['data_json'] = ''
    mod_base.insert_reaction(rnx_obj)


if __name__=='__main__':
    run_mods()