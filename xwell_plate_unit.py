import pandas as pd
from nicegui import events, ui, app
from oligoMass import molmassOligo as mmo
from datetime import datetime
import random
from collections import Counter
from lcms_chrom_data import chrom_dialog
from OligoMap_utils import oligomaps_search
from lcms_dialog_model import lcms_dialog


class click_azide():

    def __init__(self, oligos_sequence, amount_oe):
        self.seq = oligos_sequence
        self.amount = amount_oe
        self.tab = {}
        self.__rools_protocol()

    def __react_volume(self):
        o = mmo.oligoNASequence(self.seq)
        self.amount_nmol = self.amount * 1e6 / o.getExtinction()
        self.tab['amount nmol'] = round(self.amount_nmol, 2)
        self.tab['sequence'] = self.seq
        self.tab['amount oe'] = round(self.amount, 2)
        if self.amount_nmol >= 1 and self.amount_nmol <= 20:
            return 100
        elif self.amount_nmol > 20 and self.amount_nmol <= 40:
            return 200
        elif self.amount_nmol > 40 and self.amount_nmol <= 80:
            return 400
        elif self.amount_nmol > 80 and self.amount_nmol <= 600:
            return 600
        else:
            return 700

    def __rools_protocol(self):
        react_volume = self.__react_volume()
        self.tab['react volume, ul'] = react_volume
        self.tab['azide volume, ul'] = round(self.amount_nmol * 0.15, 0)
        self.tab['Cu buffer volume, ul'] = round(react_volume * 0.67, 0)
        self.tab['activator volume, ul'] = round(react_volume * 0.02, 0)
        self.tab['water volume, ul'] = round(react_volume - self.tab['azide volume, ul'] -
                                        self.tab['Cu buffer volume, ul'] - self.tab['activator volume, ul'], 0)

    def __call__(self, *args, **kwargs):
        return self.tab


class Well:
    def __init__(self, x, y, symb, num):
        self.x = x
        self.y = y
        self.symb = symb
        self.num = num
        self.color = '#f0f0f0'
        self.dye_color = '#f1a110'
        self.strftime_format = "%Y-%m-%d"
        self.oligo_data = {}


    def load_init_row(self, row):
        self.oligo_data['init_row'] = row.copy()
        self.color = self.get_support_color(self.oligo_data['init_row']['Sequence'],
                                            self.oligo_data['init_row']['Support type'])

        #print()

        self.dye_color = self.get_dye_color(self.oligo_data['init_row']['Sequence'] +
                                            self.oligo_data['init_row']['Purif type'])


    def set_init_row(self, row):
        self.oligo_data['init_row'] = row.copy()
        self.oligo_data['init_row']['Order id'] = self.oligo_data['init_row']['#']
        s = self.oligo_data['init_row']['Sequence']
        end5, end3 = self.oligo_data['init_row']["5'-end"], self.oligo_data['init_row']["3'-end"]
        if end5.lower() != 'none':
            s = f'[{end5}]{s}'
        if end3.lower() != 'none':
            s = f'{s}[{end3}]'
        self.oligo_data['init_row']['f_sequence'] = s

        self.oligo_data['init_row']['CPG, mg'] = self.get_suport_amount(self.oligo_data['init_row']['Amount, oe'])
        self.oligo_data['init_row']['Support type'] = self.get_support_type(self.oligo_data['init_row']['f_sequence'])
        self.color = self.get_support_color(self.oligo_data['init_row']['f_sequence'])
        self.dye_color = self.get_dye_color(self.oligo_data['init_row']['f_sequence'] +
                                            self.oligo_data['init_row']['Purification'])

        purif = self.get_purif_click_type(self.oligo_data['init_row']['f_sequence'])
        self.oligo_data['init_row']['Do LCMS'] = purif['lcms']
        self.oligo_data['init_row']['Do synth'] = True
        self.oligo_data['init_row']['Do cart'] = purif['cart']
        self.oligo_data['init_row']['Do hplc'] = purif['hplc']
        self.oligo_data['init_row']['Do paag'] = False
        self.oligo_data['init_row']['Do sed'] = False
        self.oligo_data['init_row']['Do click'] = purif['click']
        self.oligo_data['init_row']['Do subl'] = True

        for done in ['Done LCMS', 'Done synth', 'Done cart', 'Done hplc', 'Done paag',
                     'Done sed', 'Done click', 'Done subl']:
            self.oligo_data['init_row'][done] = False

        self.oligo_data['init_row']['Dens, oe/ml'] = 0.
        self.oligo_data['init_row']['Vol, ml'] = 0.3
        self.oligo_data['init_row']['Purity, %'] = 50.
        self.oligo_data['init_row']['Position'] = f"{self.symb}{self.num}"
        self.oligo_data['init_row']['Sequence'] = self.oligo_data['init_row']['f_sequence']
        self.oligo_data['init_row']['Purif type'] = self.oligo_data['init_row']['Purification']
        self.oligo_data['init_row']['Date'] = datetime.now().date().strftime(self.strftime_format)
        self.oligo_data['init_row']['Scale, OE'] = self.oligo_data['init_row']['Amount, oe']
        self.oligo_data['init_row']['Status'] = 'in queue'
        self.oligo_data['init_row']['DONE'] = False
        self.oligo_data['init_row']['Wasted'] = False
        self.oligo_data['init_row']['Send'] = False


    def __str__(self):
        return f"pos: {self.symb}{self.num}; coord({self.x},{self.y})"

    def get_suport_amount(self, amount):
        if amount in ['1-3', '1-5', '3-5', '5-10', '3-10']:
            return '3'
        elif amount in ['10-15']:
            return '5'
        elif amount in ['15-20', '10-20']:
            return '8'
        else:
            return '10'


    def get_support_type(self, sequence):
        oligo = mmo.oligoNASequence(sequence)
        if sequence.find('BHQ1') > -1:
            return 'bhq1_1000_hg'
        elif sequence.find('BHQ2') > -1:
            return 'bhq2_1000_hg'
        elif sequence.find('BHQ3') > -1:
            return 'bhq3_500'

        if oligo.size() <= 30:
            return 'biocomma_500'
        elif oligo.size() > 30 and oligo.size() <= 65:
            return 'biocomma_1000'
        elif oligo.size() > 65 and oligo.size() <= 110:
            return 'biocomma_2000'
        elif oligo.size() > 110 and oligo.size() <= 500:
            return 'biocomma_3000'
        else:
            return 'biocomma_3000'


    def get_support_color(self, sequence, support_type='_500'):
        if sequence.find('BHQ1') > -1:
            return '#f100a1'
        elif sequence.find('BHQ2') > -1:
            return '#a100f1'
        elif sequence.find('BHQ3') > -1:
            return '#0000f1'
        else:
        #    return '#f1f1f1'
            stype = support_type
            if stype.find('_500') > -1:
                return '#afafaf'
            if stype.find('_1000') > -1:
                return '#afaf6f'
            if stype.find('_2000') > -1:
                return '#6fafaf'
            if stype.find('_3000') > -1:
                return '#c1c1c1'
            else:
                return '#ffffff'


    def get_dye_color(self, sequence):
        if (sequence.find('FAM') > -1) or (sequence.find('6FAM') > -1):
            return '#10fa10'
        elif (sequence.find('HEX') > -1) or (sequence.find('R6G') > -1) or (sequence.find('SIMA') > -1):
            return '#f1f101'
        elif (sequence.find('VIC') > -1) or (sequence.find('JOE') > -1) or (sequence.find('Cy3') > -1):
            return '#f1f101'
        elif (sequence.find('ROX') > -1) or (sequence.find('TAMRA') > -1):
            return '#f100f1'
        elif (sequence.find('Cy5') > -1) or (sequence.find('Cy5.5') > -1) or (sequence.find('Cy7') > -1):
            return '#1a00f1'
        else:
            return '#f0f0f0'


    def get_purif_click_type(self, sequence):
        if sequence.find('FAM]') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': False}
        elif sequence.find('[HEX') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[R6G') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[SIMA') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[ROX') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[Cy') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[TAMR') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        elif sequence.find('[Alk') > -1:
            return {'cart': False, 'hplc': True, 'lcms': True, 'click': True}
        else:
            return {'cart': True, 'hplc': False, 'lcms': False, 'click': False}


def random_color_hex(r_min=0, r_max=255, g_min=0, g_max=255, b_min=0, b_max=255):
    r = random.randint(r_min, r_max)
    g = random.randint(g_min, g_max)
    b = random.randint(b_min, b_max)
    return f'#{r:02X}{g:02X}{b:02X}'


class XWell_plate:
    def __init__(self):
        self.loaded_map_id = ''
        self.plate_bkg = 'plate_bkg_2.png'
        self.well_rad = 45
        self.mouse_down = False
        self.wells = {}
        self.selected_wells = {}
        self.layer_selector = 'Base layer'
        self.order_colors = {}
        self.order_label = {}
        for x, num in zip(range(130, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
            for y, symb in zip(range(130, 900, 100), 'A B C D E F G H'.split(' ')):
                self.wells[(x, y)] = Well(x, y, symb, num)

        self.image = ui.interactive_image(f'/img/{self.plate_bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                     events=['mousedown', 'mousemove', 'mouseup'])

        with ui.context_menu() as self.context_menu:
            ui.menu_item('Редактировать хроматограмму', on_click=self.on_edit_chrom)
            ui.menu_item('Редактировать LCMS', on_click=self.on_edit_lcms)
            ui.menu_item('Edit', on_click=self.on_edit)
            #ui.menu_item('Delete', on_click=self.on_delete)

        #self.image.on('contextmenu.prevent', self.open_menu)

        self.highlight_well_select = self.image.add_layer()
        self.oligo_data_layer = self.image.add_layer()
        self.highlight_well = self.image.add_layer()
        self.draw_wells()
        self.coord_context = self.get_coord(0, 0)
        self.chrom_editor = chrom_dialog([{}])
        self.lcms_editor = lcms_dialog()
        #self.chrom_editor.on_send_chrom_data = self.on_save_chrom_data_to_base

    def open_menu(self, e):
        coord = self.get_coord(e.image_x, e.image_y)
        if coord['key'] != (0, 0):
            if 'init_row' in list(self.wells[coord['key']].oligo_data.keys()):
                self.chrom_editor.rowdata = [self.wells[coord['key']].oligo_data['init_row']]
                self.lcms_editor.rowdata = [self.wells[coord['key']].oligo_data['init_row']]
            else:
                self.chrom_editor.rowdata = []
                self.lcms_editor.rowdata = []
        else:
            self.chrom_editor.rowdata = []
            self.lcms_editor.rowdata = []
        self.context_menu.open()

    def on_edit_chrom(self):
        if self.chrom_editor.rowdata != []:
            ip = app.storage.general.get('db_IP')
            port = app.storage.general.get('db_port')
            omap = oligomaps_search(ip, port)
            omap.pincode = app.storage.user.get('pincode')
            self.chrom_editor.data_from_base = omap.search_chrom_data(int(self.chrom_editor.rowdata[0]['Order id']),
                                   self.chrom_editor.rowdata[0]['Position'])

            if self.chrom_editor.data_from_base == {}:
                self.chrom_editor.set_data_to_form_rowdata()
            else:
                self.chrom_editor.set_data_to_form_base()

            self.chrom_editor.dialog.open()

    def on_edit_lcms(self):
        if self.lcms_editor.rowdata != []:
            self.lcms_editor.set_rowdata()
            self.lcms_editor.dialog.open()

    def on_edit(self):
        ui.notify('Edit selected')

    def on_delete(self):
        ui.notify('Delete selected')


    def clear_wells(self):
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(self.wells[key].oligo_data.keys()):
                self.wells[key].oligo_data = {}
        self.selected_wells = {}
        self.draw_selected_wells()
        #self.draw_oligo_data_layer()
        self.draw_layers(self.layer_selector)


    def clear_selected_wells(self):
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(self.wells[key].oligo_data.keys()):
                if key in list(self.selected_wells.keys()):
                    self.wells[key].oligo_data = {}
        self.selected_wells = {}
        self.draw_selected_wells()
        #self.draw_oligo_data_layer()
        self.draw_layers(self.layer_selector)


    def replace_to_selected_wells(self):
        wells_list = []
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(self.wells[key].oligo_data.keys()):
                wells_list.append(well.oligo_data.copy())

        for key, well in zip(self.wells.keys(), self.wells.values()):
            self.wells[key].oligo_data = {}

        index = 0
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if index < len(wells_list):
                wells_list[index]['init_row']['Position'] = f"{well.symb}{well.num}"
                self.wells[key].load_init_row(wells_list[index]['init_row'])
                index += 1

        self.draw_selected_wells()
        self.draw_layers(self.layer_selector)

    def get_map_id(self):
        map_id = ''
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                map_id = well.oligo_data['init_row']['map #']
                if map_id != '':
                    break
        return map_id

    def well_selection_by_position(self, rowData):
        df = pd.DataFrame(rowData)
        pos_list = list(df['Position'])
        #print(pos_list)
        self.selected_wells = {}
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if f"{well.symb}{well.num}" in pos_list:
                self.selected_wells[key] = well
        return self.selected_wells


    def total_select(self):
        for key in self.wells.keys():
            self.selected_wells[key] = self.wells[key]


    def get_rowData(self):
        out = []
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                out.append(well.oligo_data['init_row'])
        return out


    def add_selrows(self, selRows, number_of_copies):
        if len(self.selected_wells.keys()) == 0:
            self.total_select()
        row_index, number = 0, 0
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if row_index < len(selRows):
                self.wells[key].set_init_row(selRows[row_index])
                number += 1
                if number == number_of_copies:
                    row_index += 1
                    number = 0
            else:
                break

        #self.draw_oligo_data_layer()
        self.draw_layers(self.layer_selector)


    def load_selrows(self, selRows):
        #print(selRows)
        self.clear_wells()
        self.well_selection_by_position(selRows)
        #self.total_select()

        row_index = 0
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if row_index < len(selRows):
                self.wells[key].load_init_row(selRows[row_index])
                row_index += 1
            else:
                break

        #self.draw_oligo_data_layer()
        self.draw_layers(self.layer_selector)


    def get_copy(self):
        return {
            'layer_selector': self.layer_selector,
            'plate_bkg': self.plate_bkg,
            'well_rad': self.well_rad,
            'mouse_down': self.mouse_down,
            'order_colors': self.order_colors,
            'order_label': self.order_label,
            'wells': self.wells.copy(),
            'selected_wells': self.selected_wells.copy(),
        }


    def set_copy(self, content):
        self.plate_bkg = content['plate_bkg']
        self.well_rad = content['well_rad']
        self.mouse_down = content['mouse_down']
        self.wells = content['wells'].copy()
        self.selected_wells = content['selected_wells'].copy()

        self.layer_selector = content['layer_selector']
        self.order_colors = content['order_colors']
        self.order_label = content['order_label']

        self.draw_selected_wells()
        #self.draw_oligo_data_layer()
        self.draw_layers(self.layer_selector)

    def get_coord(self, x, y):
        out = {'well': ('', ''), 'key': (0, 0)}
        for key in self.wells.keys():
            r2 = (x - key[0]) ** 2 + (y - key[1]) ** 2
            if r2 <= self.well_rad ** 2:
                out['well'] = self.wells[key]
                out['key'] = key
                break
        return out


    def get_selections_by_rect(self, coord1, coord2):
        if not self.ctrl_pushed:
            self.selected_wells = {}
            self.highlight_well_select.content = ""
        for key in self.wells.keys():
            if (key[0] >= coord1['key'][0]) and (key[0] <= coord2['key'][0]):
                if (key[1] >= coord1['key'][1]) and (key[1] <= coord2['key'][1]):
                    self.selected_wells[key] = self.wells[key]


    def draw_wells(self):
        for x in range(130, 1300, 100):
            for y in range(130, 900, 100):
                self.image.content += (f'<circle cx="{x}" cy="{y}" r="45" fill="gray" '
                                       f'stroke="lightgreen" stroke-width="6" />')

        for x, num in zip(range(115, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
            self.image.content += f'<text x="{x}" y="{70}" fill="orange" font-size="48">{num}</text>'
            self.image.content += f'<text x="{x}" y="{920}" fill="lightblue" font-size="48">{num}</text>'

        for y, num in zip(range(145, 900, 100), 'A B C D E F G H'.split(' ')):
            self.image.content += f'<text x="{40}" y="{y}" fill="orange" font-size="48">{num}</text>'
            self.image.content += f'<text x="{1300}" y="{y}" fill="lightblue" font-size="48">{num}</text>'


    def draw_layers(self, selector):
        if selector == 'Base layer':
            self.draw_oligo_data_layer()
        elif selector == 'Status layer':
            self.draw_status_layer()
        elif selector == 'Purification layer':
            self.draw_purification_layer()
        elif selector == 'Support layer':
            self.draw_support_layer()
        elif selector == 'Click layer':
            self.draw_click_layer()
        elif selector == 'Order layer':
            self.draw_order_layer()
        elif selector == 'Chrom layer':
            self.draw_chrom_data_layer()
        elif selector == 'LCMS layer':
            self.draw_lcms_data_layer()


    def draw_chrom_data_layer(self):
        self.oligo_data_layer.content = ""
        ip = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        omap = oligomaps_search(ip, port)
        omap.pincode = app.storage.user.get('pincode')
        ret = omap.check_chrom_data_in_base(self.loaded_map_id)
        if ret.status_code == 200:
            data = ret.json()
            for key, well in zip(self.wells.keys(), self.wells.values()):
                if 'init_row' in list(well.oligo_data.keys()):
                    if (well.oligo_data['init_row']['Order id'] in data['id'] and
                            well.oligo_data['init_row']['Position'] in data['pos']):
                        self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{well.color}" '
                                                      f'stroke="blue" stroke-width="6" />')


    def draw_lcms_data_layer(self):
        self.oligo_data_layer.content = ""
        ip = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        omap = oligomaps_search(ip, port)
        omap.pincode = app.storage.user.get('pincode')
        if self.loaded_map_id != '':
            pos_list = omap.check_lcms_data_in_base(int(self.loaded_map_id))
        else:
            pos_list = []
        if len(pos_list) > 0:
            for key, well in zip(self.wells.keys(), self.wells.values()):
                if 'init_row' in list(well.oligo_data.keys()):
                    if well.oligo_data['init_row']['Position'] in pos_list:
                        self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{well.color}" '
                                                      f'stroke="blue" stroke-width="6" />')


    def draw_oligo_data_layer(self):
        self.oligo_data_layer.content = ""
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{well.color}" '
                                                  f'stroke="{well.color}" stroke-width="6" />')
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="20" fill="{well.dye_color}" '
                                                  f'stroke="{well.dye_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Order id']} ({well.symb}{well.num})"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 25}" y="{well.y - 15}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')
                if 'Support type' in list(well.oligo_data['init_row'].keys()):
                    support_type = well.oligo_data['init_row']['Support type']
                    support = well.oligo_data['init_row']['Support type'][support_type.find('_'):]
                else:
                    support = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 35}" y="{well.y - 3}" '
                                                  f'fill="black" font-size="14">size{support}</text>')

    def draw_support_layer(self):
        self.oligo_data_layer.content = ""
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{well.color}" '
                                                  f'stroke="{well.color}" stroke-width="6" />')
                #self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="20" fill="{well.dye_color}" '
                #                                  f'stroke="{well.dye_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Order id']}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 25}" y="{well.y - 15}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')
                if 'Support type' in list(well.oligo_data['init_row'].keys()):
                    support_type = well.oligo_data['init_row']['Support type']
                    if support_type.find('bhq') == -1:
                        support = well.oligo_data['init_row']['Support type'][support_type.find('_'):]
                        support = f'size{support}'
                    else:
                        support = well.oligo_data['init_row']['Support type'][support_type.find('bhq'):4]
                        support = support.upper()
                else:
                    support = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 35}" y="{well.y + 5}" '
                                                  f'fill="black" font-size="14">{support}</text>')
                if 'CPG, mg' in list(well.oligo_data['init_row'].keys()):
                    amount = well.oligo_data['init_row']['CPG, mg']
                else:
                    amount = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 20}" y="{well.y + 20}" '
                                                  f'fill="black" font-size="14">{amount}</text>')


    def set_repeated_wells_color(self):
        self.repeated_wells_colors = {}
        repeated_orders = []#Counter([well.oligo_data['init_row']['Order id'] for well in self.wells.values()])
        color = '#1f1fff'
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                repeated_orders.append(well.oligo_data['init_row']['Order id'])
                self.repeated_wells_colors[well.oligo_data['init_row']['Order id']] = '#000000'
        repeated_orders = Counter(repeated_orders)
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                if repeated_orders[well.oligo_data['init_row']['Order id']] > 1:
                    #color = random_color_hex(150, 250, 150, 250, 150, 250)
                    if self.repeated_wells_colors[well.oligo_data['init_row']['Order id']] == '#000000':
                        if color == '#1f1fff':
                            color = 'orange'
                        else:
                            color = '#1f1fff'
                        self.repeated_wells_colors[well.oligo_data['init_row']['Order id']] = color

    def draw_purification_layer(self):
        self.oligo_data_layer.content = ""
        self.set_repeated_wells_color()
        for key, well in zip(self.wells.keys(), self.wells.values()):
            color = self.get_purification_color(well)
            pur_type = self.get_purification_text(well)
            if 'init_row' in list(well.oligo_data.keys()):
                well_color = self.repeated_wells_colors[well.oligo_data['init_row']['Order id']]
                if well_color == '#000000':
                    well_color = color
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="30" fill="{color}" '
                                                  f'stroke="{color}" stroke-width="6" />')
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="15" fill="{well_color}" '
                                                  f'stroke="{well_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Order id']}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 25}" y="{well.y - 15}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')

                self.oligo_data_layer.content += (f'<text x="{well.x - 30}" y="{well.y + 5}" '
                                                  f'fill="black" font-size="14">{pur_type}</text>')


    def set_invoce_data(self, order_tab):
        data = pd.DataFrame(order_tab)
        self.order_colors = {}
        self.order_label = {}
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                self.wells[key].oligo_data['invoce'] = (
                    data[data['#'] == well.oligo_data['init_row']['Order id']]['order id'].max())
                self.order_colors[self.wells[key].oligo_data['invoce']] = '#f1f1f1'
                text_key = self.wells[key].oligo_data['invoce']
                if text_key.find('УТ') > -1:
                    self.order_label[self.wells[key].oligo_data['invoce']] = text_key[text_key.find('УТ'):]
                else:
                    self.order_label[self.wells[key].oligo_data['invoce']] = text_key[:5]

        for order_id in list(self.order_colors.keys()):
            color = random_color_hex(100, 250, 100, 250, 100, 250)
            self.order_colors[order_id] = color

    def draw_order_layer(self):
        self.oligo_data_layer.content = ""
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                if 'invoce' in list(well.oligo_data.keys()):
                    color = self.order_colors[well.oligo_data['invoce']]
                else:
                    color = 'red'
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{color}" '
                                                  f'stroke="{color}" stroke-width="6" />')
                #self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="20" fill="{well.dye_color}" '
                #                                  f'stroke="{well.dye_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Order id']}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 20}" y="{well.y - 15}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')

                if 'invoce' in list(well.oligo_data.keys()):
                    show_data = f"{self.order_label[well.oligo_data['invoce']]}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 33}" y="{well.y + 4}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')



    def get_status_color(self, status):
        if status == 'finished':
            return 'lightgreen'
        elif status == 'in queue' or status == 'synthesis':
            return 'orange'
        else:
            return 'white'

    def get_purification_color(self, well):
        color = 'white'
        if 'init_row' in list(well.oligo_data.keys()):
            if well.oligo_data['init_row']['Do cart']:
                color = 'lightgreen'
            else:
                color = '#a100f1'
            if well.oligo_data['init_row']['Wasted']:
                color = 'red'
        return color

    def get_status_color_well(self, well):
        color = 'lightgreen'
        if 'init_row' in list(well.oligo_data.keys()):
            if well.oligo_data['init_row']['Status'] != 'finished':
                color = 'orange'
            if well.oligo_data['init_row']['Wasted']:
                color = 'red'
        return color

    def get_purification_text(self, well):
        ret = 'RP-Cart'
        if 'init_row' in list(well.oligo_data.keys()):
            if well.oligo_data['init_row']['Do cart']:
                ret = 'RP-Cart'
            elif well.oligo_data['init_row']['Do sed']:
                ret = 'Sediment'
            else:
                ret = 'HPLC'
            if well.oligo_data['init_row']['Wasted']:
                ret = 'Wasted'
            if well.oligo_data['init_row']['Do paag']:
                ret = 'PAAG'
        return ret

    def draw_status_layer(self):
        self.oligo_data_layer.content = ""

        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                color = self.get_status_color_well(well)
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{color}" '
                                                  f'stroke="{color}" stroke-width="6" />')
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="10" fill="{well.dye_color}" '
                                                  f'stroke="{well.dye_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Status' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Status']}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 30}" y="{well.y}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')
                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Order id']}"
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 20}" y="{well.y - 15}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')


    def culc_click(self):
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                seq = well.oligo_data['init_row']['Sequence']
                if seq.find('[Alk]') > -1:
                    amount = well.oligo_data['init_row']['Dens, oe/ml'] * well.oligo_data['init_row']['Vol, ml']
                    self.wells[key].oligo_data['click_data'] = click_azide(seq, amount)()
                    self.wells[key].oligo_data['click_data']['water volume, ul'] = round(
                    self.wells[key].oligo_data['click_data']['water volume, ul'] -
                    well.oligo_data['init_row']['Vol, ml'] * 1000, 0)


    def draw_click_layer(self):
        self.oligo_data_layer.content = ""
        self.culc_click()
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                color = self.get_purification_color(well)
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="32" fill="{well.color}" '
                                                  f'stroke="{well.color}" stroke-width="6" />')
                self.oligo_data_layer.content += (f'<circle cx="{well.x}" cy="{well.y}" r="10" fill="{well.dye_color}" '
                                                  f'stroke="{well.dye_color}" stroke-width="6" />')

                #self.oligo_data_layer.content += f'<path d = "M {well.x},{well.y} a 80,80 0 0,1 160,0" fill = "{well.dye_color}" / >'

                if 'Purif type' in list(well.oligo_data['init_row'].keys()):
                    show_data = f"{well.oligo_data['init_row']['Purif type']}"
                    if show_data.find('_') > -1:
                        show_data = show_data[show_data.find('_') + 1:]
                    else:
                        show_data = ''
                else:
                    show_data = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 15}" y="{well.y - 18}" '
                                                  f'fill="black" font-size="16">{show_data}</text>')
                if 'click_data' in list(well.oligo_data.keys()):
                    if 'azide volume, ul' in list(well.oligo_data['click_data'].keys()):
                        show_data = f"az{int(well.oligo_data['click_data']['azide volume, ul'])}"
                    else:
                        show_data = ''
                    self.oligo_data_layer.content += (f'<text x="{well.x - 35}" y="{well.y + 0}" '
                                                      f'fill="black" font-size="16">{show_data}</text>')

                    if 'activator volume, ul' in list(well.oligo_data['click_data'].keys()):
                        show_data = f"ac{int(well.oligo_data['click_data']['activator volume, ul'])}"
                    else:
                        show_data = ''
                    self.oligo_data_layer.content += (f'<text x="{well.x + 5}" y="{well.y + 0}" '
                                                      f'fill="black" font-size="16">{show_data}</text>')

                    if 'Cu buffer volume, ul' in list(well.oligo_data['click_data'].keys()):
                        show_data = f"b{int(well.oligo_data['click_data']['Cu buffer volume, ul'])}"
                    else:
                        show_data = ''
                    self.oligo_data_layer.content += (f'<text x="{well.x - 35}" y="{well.y + 23}" '
                                                      f'fill="black" font-size="16">{show_data}</text>')

                    if 'water volume, ul' in list(well.oligo_data['click_data'].keys()):
                        show_data = f"w{int(well.oligo_data['click_data']['water volume, ul'])}"
                    else:
                        show_data = ''
                    self.oligo_data_layer.content += (f'<text x="{well.x + 5}" y="{well.y + 23}" '
                                                      f'fill="black" font-size="16">{show_data}</text>')

                #self.oligo_data_layer.content += (f'<text x="{well.x - 20}" y="{well.y}" '
                #                                  f'fill="black" font-size="16">{show_data}</text>')


    def draw_selected_wells(self):
        self.highlight_well_select.content = ""
        for key in self.selected_wells.keys():
            if self.selected_wells[key] != ('', ''):
                x = self.selected_wells[key].x
                y = self.selected_wells[key].y
                self.highlight_well_select.content += (f"<circle cx='{x}' cy='{y}' "
                                    f"r='45' fill='yellow' stroke='lightgreen' stroke-width='6' />")


    def get_selected_pos_list(self):
        pos_list = []
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            pos_list.append(f"{well.symb}{well.num}")
        return pos_list

    def get_selected_rows(self):
        lb_list = []
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                lb_list.append(well.oligo_data['init_row'])
        return lb_list

    def get_selected_labels(self):
        lb_list = []
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if 'init_row' in list(well.oligo_data.keys()):
                lb_list.append(well.oligo_data['init_row'])
        pass_df = self.print_pass(lb_list)
        return pass_df


    def  print_pass(self, rowData):
        out_tab = []
        index_ = 1

        for row in rowData:
            nseq = row['Sequence']
            if type(row['Purif type']) == str:
                if row['Purif type'].find('_') > 0:
                    fluoro = row['Purif type'][row['Purif type'].find('_')+1:]
                    nseq = row['Sequence'].replace('[Alk]', f"[{fluoro}]")
                    #print(nseq)
                    o = mmo.oligoNASequence(nseq)
                else:
                    o = mmo.oligoNASequence(row['Sequence'])
            d = {}
            d['#'] = index_
            index_ += 1
            #d['Position'] = row['Position']
            d['Name'] = row['Name'] + f"  ({row['Synt number']}_{row['Position']})"
            d['Sequence'] = nseq
            d['Amount,_oe'] = int(round(float(row['Dens, oe/ml']) * float(row['Vol, ml']), 0))
            if o.getExtinction() > 0:
                d['Amount,_nmol'] = int(round(d['Amount,_oe'] * 1e6 / o.getExtinction(), 0))
            else:
                d['Amount,_nmol'] = 0.
            d['Desolving'] = int(d['Amount,_nmol'] * 10)

            d['Purification'] = row['Purif type']
            d['order_ID'] = row['Order id']
            d['Status'] = row['Status']
            try:
                d['Mass,_Da'] = round(o.getAvgMass(), 2)
            except:
                d['Mass,_Da'] = 'unknown modiff'
            d['Extinction'] = o.getExtinction()
            #print(d)
            out_tab.append(d)
        return pd.DataFrame(out_tab)


    def mouse_handler(self, e: events.MouseEventArguments):
        if e.button == 0:
            self.mouse_handler_0(e)
        elif e.button == 2:
            self.open_menu(e)

    def mouse_handler_0(self, e: events.MouseEventArguments):

        ip = app.storage.user.get('client_ip')

        coord = self.get_coord(e.image_x, e.image_y)

        self.ctrl_pushed = e.ctrl

        if e.type == 'mousedown':
            self.mouse_down = True
            self.init_coord = coord.copy()
            if not self.ctrl_pushed:
                self.selected_wells = {}
                self.highlight_well_select.content = ""


        if e.type == 'mousemove':
            if coord['well'] != ('', '') and not self.mouse_down:
                self.highlight_well.content = (f"<circle cx='{coord['key'][0]}' cy='{coord['key'][1]}' "
                                          f"r='45' fill='lightblue' stroke='lightgreen' stroke-width='6' />")
            if self.mouse_down:
                if coord['well'] != ('', ''):
                    self.move_coord = coord.copy()
                    #self.selected_wells[coord['key']] = coord['well']
                    self.get_selections_by_rect(self.init_coord, self.move_coord)
                    self.draw_selected_wells()

        if e.type == 'mouseup':
            self.mouse_down = False
            if self.ctrl_pushed:
                if coord['key'] in list(self.selected_wells.keys()):
                    self.selected_wells.pop(coord['key'])
                    self.highlight_well_select.content = ""
                    self.draw_selected_wells()
                else:
                    self.selected_wells[coord['key']] = coord['well']
                    self.draw_selected_wells()
