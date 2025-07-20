from nicegui import events, ui, app

class Well:
    def __init__(self, x, y, symb, num):
        self.x = x
        self.y = y
        self.symb = symb
        self.num = num

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
        self.highlight_well = self.image.add_layer()
        self.highlight_well_select = self.image.add_layer()
        self.draw_wells()


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
                self.image.content += f'<circle cx="{x}" cy="{y}" r="45" fill="gray" stroke="lightgreen" stroke-width="6" />'

        for x, num in zip(range(115, 1300, 100), '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')):
            self.image.content += f'<text x="{x}" y="{70}" fill="orange" font-size="48">{num}</text>'
            self.image.content += f'<text x="{x}" y="{920}" fill="lightblue" font-size="48">{num}</text>'

        for y, num in zip(range(145, 900, 100), 'A B C D E F G H'.split(' ')):
            self.image.content += f'<text x="{40}" y="{y}" fill="orange" font-size="48">{num}</text>'
            self.image.content += f'<text x="{1300}" y="{y}" fill="lightblue" font-size="48">{num}</text>'


    def draw_selected_wells(self):
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
