from nicegui import events, ui, app
from oligoMass import molmassOligo as mmo
from datetime import datetime

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
        self.color = self.get_support_color(self.oligo_data['init_row']['Sequence'])
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


    def get_support_color(self, sequence):
        if sequence.find('BHQ1') > -1:
            return '#f100a1'
        elif sequence.find('BHQ2') > -1:
            return '#a100f1'
        elif sequence.find('BHQ3') > -1:
            return '#0000f1'
        else:
        #    return '#f1f1f1'
            stype = self.get_support_type(sequence)
            if stype.find('_500') > -1:
                return '#313131'
            if stype.find('_1000') > -1:
                return '#717171'
            if stype.find('_2000') > -1:
                return '#919191'
            if stype.find('_3000') > -1:
                return '#a1a1a1'
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
        elif (sequence.find('Cy5') > -1) or (sequence.find('Cy5.5') > -1):
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



class XWell_plate:
    def __init__(self, frontend, backend):
        self.frontend = frontend
        self.backend = backend
        self.plate_bkg = 'plate_bkg_2.png'
        self.well_rad = 45
        self.mouse_down = False
        self.wells = {}
        self.selected_wells = {}
        for x, num in zip(range(130, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
            for y, symb in zip(range(130, 900, 100), 'A B C D E F G H'.split(' ')):
                self.wells[(x, y)] = Well(x, y, symb, num)

        self.image = ui.interactive_image(f'/img/{self.plate_bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                     events=['mousedown', 'mousemove', 'mouseup'])
        self.highlight_well_select = self.image.add_layer()
        self.oligo_data_layer = self.image.add_layer()
        self.highlight_well = self.image.add_layer()
        self.draw_wells()


    def clear_wells(self):
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(self.wells[key].oligo_data.keys()):
                self.wells[key].oligo_data = {}
        self.selected_wells = {}
        self.draw_selected_wells()
        self.draw_oligo_data_layer()


    def clear_selected_wells(self):
        for key, well in zip(self.wells.keys(), self.wells.values()):
            if 'init_row' in list(self.wells[key].oligo_data.keys()):
                if key in list(self.selected_wells.keys()):
                    self.wells[key].oligo_data = {}
        self.selected_wells = {}
        self.draw_selected_wells()
        self.draw_oligo_data_layer()


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

        self.draw_oligo_data_layer()


    def load_selrows(self, selRows):
        self.clear_wells()
        self.total_select()

        row_index = 0
        for key, well in zip(self.selected_wells.keys(), self.selected_wells.values()):
            if row_index < len(selRows):
                self.wells[key].load_init_row(selRows[row_index])
                row_index += 1
            else:
                break

        self.draw_oligo_data_layer()


    def get_copy(self):
        return {
            'plate_bkg': self.plate_bkg,
            'well_rad': self.well_rad,
            'mouse_down': self.mouse_down,
            'wells': self.wells.copy(),
            'selected_wells': self.selected_wells.copy(),
        }


    def set_copy(self, content):
        self.plate_bkg = content['plate_bkg']
        self.well_rad = content['well_rad']
        self.mouse_down = content['mouse_down']
        self.wells = content['wells'].copy()
        self.selected_wells = content['selected_wells'].copy()

        self.draw_selected_wells()
        self.draw_oligo_data_layer()


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
                if 'Order id' in list(well.oligo_data['init_row'].keys()):
                    support_type = well.oligo_data['init_row']['Support type']
                    support = well.oligo_data['init_row']['Support type'][support_type.find('_'):]
                else:
                    support = ''
                self.oligo_data_layer.content += (f'<text x="{well.x - 35}" y="{well.y - 3}" '
                                                  f'fill="black" font-size="14">size{support}</text>')


    def draw_selected_wells(self):
        self.highlight_well_select.content = ""
        for key in self.selected_wells.keys():
            if self.selected_wells[key] != ('', ''):
                x = self.selected_wells[key].x
                y = self.selected_wells[key].y
                self.highlight_well_select.content += (f"<circle cx='{x}' cy='{y}' "
                                    f"r='45' fill='yellow' stroke='lightgreen' stroke-width='6' />")


    def mouse_handler(self, e: events.MouseEventArguments):

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

            self.backend.client_frontend[ip] = self.frontend.get_model()
            #print('GET MODEL', self.backend.client_frontend[ip])
