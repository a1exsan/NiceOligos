from io import BytesIO
from PIL import Image, ImageDraw
from nicegui import ui, app
import base64

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

import hashlib
from datetime import datetime

import random
import math


class moleculeInfo():
    def __init__(self, text_struct):
        self.text_struct = text_struct
        if text_struct.find('InChI') > -1:
            self.mol = Chem.inchi.MolFromInchi(text_struct)
        else:
            self.mol = Chem.MolFromSmiles(text_struct)
        if self.mol is None:
            raise ValueError('Некорректный InChI')

    def draw_svg(self, width=600, height=400):
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        rdMolDraw2D.SetDarkMode(drawer)
        drawer.DrawMolecule(self.mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def get_props(self):
        return {
            'Mol weight, Da ': round(Descriptors.MolWt(self.mol), 2),
            'Brutto ': rdMolDescriptors.CalcMolFormula(self.mol)
        }



class image_background():
    def __init__(self, width, height, color='black'):
        self.width = width
        self.height = height
        self.color = color

        self.background = self.generate_background()
        buffer = BytesIO()
        self.background.save(buffer, format='PNG')
        buffer.seek(0)
        self.background_base64 = base64.b64encode(buffer.read()).decode()

    def generate_background(self):
        img = Image.new('RGB', (self.width, self.height), color=self.color)
        draw = ImageDraw.Draw(img)
        draw.rectangle([(0, 0), (self.width, self.height)], fill=self.color)
        return img


class molecule():
    def __init__(self, id, struct_data, image):
        self.id = id
        self.struct_data = struct_data
        self.struct = moleculeInfo(struct_data)

        self.image = image
        self.layer = self.image.add_layer()
        self.width = 300
        self.height = 300
        self.x_pos = 300
        self.y_pos = 300
        self.x_center = self.x_pos + self.width // 2
        self.y_center = self.y_pos + self.height // 2
        self.radius = self.width // 3
        self.is_hold = False

        self.svg = self.struct.draw_svg(width=self.width, height=self.height)

    def set_coord(self, x, y):
        self.x_pos = x
        self.y_pos = y
        self.x_center = self.x_pos + self.width // 2
        self.y_center = self.y_pos + self.height // 2

    def isin(self, x, y):
        if ((x - self.x_center)**2 + (y - self.y_center)**2) <= self.radius**2:
            return True
        else:
            return False

    def on_hold(self, e):
        x, y = e.image_x, e.image_y
        if self.isin(x, y):
            self.is_hold = True
            self.hold_x, self.hold_y = x, y

    def on_move(self, e):
        x, y = e.image_x, e.image_y
        if self.is_hold:
            self.x_pos += x - self.hold_x
            self.y_pos += y - self.hold_y
            self.x_center = self.x_pos + self.width // 2
            self.y_center = self.y_pos + self.height // 2
            self.hold_x, self.hold_y = x, y
            self.draw()


    def on_up(self, e):
        x, y = e.image_x, e.image_y
        self.is_hold = False


    def draw(self):
        fillop = 0.2
        self.layer.content = ""
        self.layer.content += f'<g transform="translate({self.x_pos}, {self.y_pos})">{self.svg}</g>'
        self.layer.content += (f'<circle cx="{self.x_center}" cy="{self.y_center}" fill-opacity="{fillop}"'
                                    f'r="{self.radius}" fill="lightgreen" stroke="gray" stroke-width="0" />')



class molecule_interaction():
    def __init__(self, molecules):
        self.molecules = molecules
        self.inter_trehold = 0.2

    def gen_interaction_matrix(self):
        inter_dict = {}
        for i, mol1 in enumerate(self.molecules.values()):
            for j, mol2 in enumerate(self.molecules.values()):
                if i < j:
                    d = ((mol1.x_center - mol2.x_center)**2 + (mol1.y_center - mol2.y_center)**2)**0.5
                    r1 = mol1.radius
                    r2 = mol2.radius
                    inter_area = self.intersection_area(r1, r2, d)
                    if inter_area / (3.14 * r1**2) >= self.inter_trehold:
                        inter_dict[(mol1.id, mol2.id)] = {'inter_area': inter_area}
        return inter_dict


    def intersection_area(self, R1, R2, d):
        # Если окружности не пересекаются
        if d >= R1 + R2:
            return 0.0
        # Если одна окружность полностью внутри другой
        if d <= abs(R1 - R2):
            return math.pi * min(R1, R2) ** 2

        # Вычисляем площадь пересечения
        part1 = R1 ** 2 * math.acos((d ** 2 + R1 ** 2 - R2 ** 2) / (2 * d * R1))
        part2 = R2 ** 2 * math.acos((d ** 2 + R2 ** 2 - R1 ** 2) / (2 * d * R2))
        part3 = 0.5 * math.sqrt((-d + R1 + R2) * (d + R1 - R2) * (d - R1 + R2) * (d + R1 + R2))

        return part1 + part2 - part3

    def on_inter_detect(self, e):
        inter_mx = self.gen_interaction_matrix()
        print(inter_mx)


class structure_field():
    def __init__(self):
        self.start_x_pos = 300
        self.start_y_pos = 300
        self.width = 1900
        self.height = 1200
        self.bkg_color = '#2c2c2c'
        bkg = image_background(self.width, self.height, color=self.bkg_color)
        self.image = ui.interactive_image(f"data:image/png;base64,{bkg.background_base64}",
                                          size=(self.width, self.height), on_mouse=self.mouse_handler,
                                          events=['mousedown', 'mousemove', 'mouseup', 'click']
                                          )
        self.molucules = {}
        self.focus_mol = {}
        self.interaction_mx = molecule_interaction(self.molucules)


    def gen_id(self):
        now = datetime.now()
        datetime_str = now.isoformat()
        hash_object = hashlib.sha256(datetime_str.encode('utf-8'))
        hash_hex = hash_object.hexdigest()
        return hash_hex

    def draw(self):
        for mol in self.molucules.values():
            mol.draw()

    def adjust_structure(self, text_data):
        id = self.gen_id()
        self.molucules[id] = molecule(id, text_data, self.image)
        x = self.start_x_pos + random.randint(30, 60)
        y = self.start_y_pos + random.randint(30, 60)
        self.molucules[id].set_coord(x, y)
        self.start_x_pos, self.start_y_pos = x, y
        self.focus_mol = id
        self.draw()
        self.interaction_mx = molecule_interaction(self.molucules)

    def mouse_handler(self, e):
        x, y = e.image_x, e.image_y

        if e.type == 'mousedown':
            for mol in self.molucules.values():
                mol.on_hold(e)

        if e.type == 'mousemove':
            for mol in self.molucules.values():
                mol.on_move(e)
                self.interaction_mx.on_inter_detect(e)

        if e.type == 'mouseup':
            for mol in self.molucules.values():
                mol.on_up(e)

        if e.type == 'click':
            pass


def on_add_structure():
    field.adjust_structure(structure_data.value)

mol1 = "InChI=1S/C30H32N6O4/c1-5-32-24-15-26-22(12-17(24)3)28(23-13-18(4)25(33-6-2)16-27(23)40-26)20-9-8-19(14-21(20)30(38)39)29(37)34-10-7-11-35-36-31/h8-9,12-16,32H,5-7,10-11H2,1-4H3,(H,34,37)(H,38,39)/b33-25-"

with ui.row():
    field = structure_field()
    with ui.column():
        mol_name = ui.input(label='Molecule name', value='molecule_1').classes('w-[400px]')
        with ui.row():
            structure_data = ui.textarea(label='Inchi or Smiles',
                                         value=mol1,
                                         ).classes('w-[600px]')
            ui.button('Add structure', on_click=on_add_structure)

ui.run(port=8081)
