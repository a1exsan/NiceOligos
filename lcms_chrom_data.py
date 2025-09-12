from nicegui import ui, events
from io import BytesIO
from PIL import Image, ImageDraw
import base64
import pandas as pd
from invoce_chart import PolylineSVG
from scipy.interpolate import interp1d
import numpy as np


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

class plot_gui_base():
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.bkg_color = '#2c2c2c'
        self.rect_fill_color = 'orange'
        self.rect_fill_opacity = '0.1'
        self.stroke_width = 2
        self.rect = {'x': 0, 'y': 0, 'w': 0, 'h': 0, 'x_down': 0, 'y_down': 0}
        self.down_btn = False
        self.hiding_rect_after_up = False
        self.draw_space = {'p1': (2, 2), 'p2': (self.width-2, self.height-2)}

        bkg = image_background(self.width, self.height, color=self.bkg_color)
        self.image = ui.interactive_image(f"data:image/png;base64,{bkg.background_base64}",
                                          size=(self.width, self.height), on_mouse=self.mouse_handler,
                                          events=['mousedown', 'mousemove', 'mouseup', 'click']
                                          )
        self.rect_layer = self.image.add_layer()


    def reset_bkg_color(self, color):
        self.bkg_color = color
        bkg = image_background(self.width, self.height, color=self.bkg_color)
        self.image.set_source(f"data:image/png;base64,{bkg.background_base64}")

    def check_mouse_on_image_area(self, x, y):
        if x <= self.draw_space['p1'][0] or x >= self.draw_space['p2'][0]:
            return False
        if y <= self.draw_space['p1'][1] or y >= self.draw_space['p2'][1]:
            return False
        return True

    def draw_rect(self):
        self.rect_layer.content = ''
        if self.down_btn:
            self.rect_layer.content += (f'<rect x={self.rect["x"]} y={self.rect["y"]} width={self.rect["w"]} '
                                        f'height={self.rect["h"]} '
                                        f' rx=5 ry=5 fill="{self.rect_fill_color}" '
                                        f'fill-opacity="{self.rect_fill_opacity}"'
                                        f'stroke="{self.rect_fill_color}" stroke-width="{self.stroke_width}"/>')


    def on_mouse_down(self, e):
        #print(e)
        pass

    def on_mouse_up(self, e):
        pass

    def on_mouse_move(self, e):
        pass

    def on_mouse_click(self, e):
        pass

    def mouse_handler(self, e):
        x, y = e.image_x, e.image_y

        if e.type == 'mousedown':
            self.down_btn = True
            self.rect['x_down'] = x
            self.rect['y_down'] = y
            self.on_mouse_down(e)

        if e.type == 'mouseup':
            self.down_btn = False
            if self.hiding_rect_after_up:
                self.draw_rect()
            self.on_mouse_up(e)

        if e.type == 'mousemove':
            if not self.check_mouse_on_image_area(x, y):
                #self.down_btn = False
                if self.hiding_rect_after_up:
                    self.draw_rect()
            if self.down_btn and self.check_mouse_on_image_area(x, y):
                self.rect['w'] = x - self.rect['x_down']
                if self.rect['w'] < 0:
                    self.rect['x'] = x
                    self.rect['w'] = abs(self.rect['w'])
                else:
                    self.rect['x'] = self.rect['x_down']

                self.rect['h'] = y - self.rect['y_down']
                if self.rect['h'] < 0:
                    self.rect['y'] = y
                    self.rect['h'] = abs(self.rect['h'])
                else:
                    self.rect['y'] = self.rect['y_down']
                self.draw_rect()
            self.on_mouse_move(e)

        if e.type == 'click':
            self.on_mouse_click(e)
            
            
class diagram_base(plot_gui_base):
    def __init__(self, width, height):
        super().__init__(width, height)
        self.x_axis_label = 'x coordinate'
        self.y_axis_label = 'y coordinate'
        self.x_axis_visible = True
        self.y_axis_visible = True
        self.axis_layer = self.image.add_layer()
        self.axis_font_size = 20
        self.axis_font_color = 'white'
        self.axis_line_color = 'white'
        self.axis_line_width = 2
        self.set_draw_space(50 + self.axis_font_size, 30, self.width - 30, self.height - 50 - self.axis_font_size)

    def set_draw_space(self, x1, y1, x2, y2):
        self.axis_lines = {}
        self.draw_space['p1'] = (x1, y1)
        self.draw_space['p2'] = (x2, y2)
        self.axis_lines['vert_x1'] = x1
        self.axis_lines['vert_y1'] = y1
        self.axis_lines['vert_x2'] = x1
        self.axis_lines['vert_y2'] = y2
        self.axis_lines['hor_x1'] = x1
        self.axis_lines['hor_y1'] = y2
        self.axis_lines['hor_x2'] = x2
        self.axis_lines['hor_y2'] = y2


    def draw_axis(self):
        self.axis_layer.content = ''

        if self.y_axis_visible:
            self.axis_layer.content += (f'<text x="{self.axis_font_size + 5}" y="{self.height // 2 + 20}" '
                                        f'transform="rotate(-90 {self.axis_font_size + 5}, {self.height // 2 + 20})"'
                                        f'fill="{self.axis_font_color}" '
                                        f'font-size="{self.axis_font_size}">{self.y_axis_label}</text>')

            self.axis_layer.content += (f'<line x1="{self.axis_lines["vert_x1"]}" y1="{self.axis_lines["vert_y1"]}" '
                                        f'x2="{self.axis_lines["vert_x2"]}" y2="{self.axis_lines["vert_y2"]}" '
                                        f'stroke="{self.axis_line_color}" '
                                        f'stroke-width="{self.axis_line_width}"/>')

        if self.x_axis_visible:
            self.axis_layer.content += (f'<text x="{self.width // 2 - 40}" '
                                        f'y="{self.height - self.axis_font_size}" '
                                        f'transform="rotate(0 {self.axis_font_size + 5}, {self.height // 2 + 20})"'
                                        f'fill="{self.axis_font_color}" '
                                        f'font-size="{self.axis_font_size}">{self.x_axis_label}</text>')

            self.axis_layer.content += (f'<line x1="{self.axis_lines["hor_x1"]}" y1="{self.axis_lines["hor_y1"]}" '
                                        f'x2="{self.axis_lines["hor_x2"]}" y2="{self.axis_lines["hor_y2"]}" '
                                        f'stroke="{self.axis_line_color}" '
                                        f'stroke-width="{self.axis_line_width}"/>')


class asc_reader():
    def __init__(self, filename):
        self.filename = filename
        if filename != '':
            self.read_file(filename)
        else:
            self.series = {}
            self.series_index = {}
            self.series_units = {}
            self.series_data = {}

    def parse(self, data):
        self.series = {}
        self.series_index = {}
        self.series_units = {}
        self.series_data = {}
        for i, row in enumerate(data):
            if i == 1:
                k_key = 'none'
                for j, key in enumerate(row.split('\t')):
                        if key not in ['']:
                            k_key = key
                            self.series_index[key] = [j]
                        else:
                            self.series_index[k_key].append(j)
            if i == 2:
                for j, key in enumerate(row.split('\t')):
                    self.series_units[j] = key
                    self.series_data[j] = []
            if i > 2:
                for j, key in enumerate(row.split('\t')):
                    if key not in ['', '\n', '\r\n']:
                        self.series_data[j].append(key)
        for key in self.series_index.keys():
            if key not in ['\n', '\r\n']:
                i1 = self.series_index[key][0]
                i2 = self.series_index[key][1]
                self.series[key] = {}
                self.series[key][self.series_units[i1]] = self.series_data[i1]
                self.series[key][self.series_units[i2]] = self.series_data[i2]


    def read_file(self, filename):
        with open(filename, 'r') as file:
            data = file.readlines()
            self.parse(data)

    def get_series_by_contains(self, contains_list):
        out_series = {}
        for key in self.series.keys():
            for label in contains_list:
                if str(key).find(label) > -1:
                    out_series[key] = self.series[key].copy()
        return out_series


class coord_scaler():
    def __init__(self, series, plot_space):
        self.plot_space = plot_space
        self.data = pd.DataFrame(series)
        keys = list(self.data.keys())
        self.points = []
        for key in self.data.keys():
            try:
                self.data[key] = self.data[key].astype(float)
            except:
                pass
        try:
            self.points_x, self.points_y = self.set_draw_points(self.data[keys[0]], self.data[keys[1]])
            self.points = [(x, y) for x, y in zip(self.points_x, self.points_y)]
        except:
            self.points_x, self.points_y = self.set_draw_points(self.data[keys[0]],
                                                                pd.Series([100 for i in range(self.data.shape[0])]))
            self.points = [(x, y) for x, y in zip(self.points_x, self.points_y)]


    def set_draw_points(self, x, y):
        x_min, y_min, x_max, y_max = min(x), min(y), max(x), max(y)

        if x_max - x_min == 0:
            x_max = (x_min + 1) * 10

        if y_max - y_min == 0:
            y_max = (y_min + 1) * 10

        x_cf = (self.plot_space['p2'][0] - self.plot_space['p1'][0]) / (x_max - x_min)
        y_cf = (self.plot_space['p2'][1] - self.plot_space['p1'][1] - 30) / (y_max - y_min)

        xx = self.plot_space['p1'][0] + x_cf * (x - x_min)
        yy = y_cf * (y - y_min)
        yy = self.plot_space['p2'][1] - yy
        return xx, yy

def set_draw_points(x, y, ref_x, ref_y, plot_space):
        x_min, y_min, x_max, y_max = min(ref_x), min(ref_y), max(ref_x), max(ref_y)

        if x_max - x_min == 0:
            x_max = (x_min + 1) * 10

        if y_max - y_min == 0:
            y_max = (y_min + 1) * 10

        x_cf = (plot_space['p2'][0] - plot_space['p1'][0]) / (x_max - x_min)
        y_cf = (plot_space['p2'][1] - plot_space['p1'][1] - 30) / (y_max - y_min)

        xx = plot_space['p1'][0] + x_cf * (x - x_min)
        yy = y_cf * (y - y_min)
        yy = plot_space['p2'][1] - yy
        return xx, yy

class interpolate_crom_line():
    def __init__(self, X, Y, point_numbers=120):
        self.points_number = point_numbers
        self.X = np.array(X)
        self.Y = np.array(Y)

    def __call__(self, *args, **kwargs):
        linear_interp = interp1d(self.X, self.Y, kind='linear')
        self.new_X = np.linspace(min(self.X), max(self.X), self.points_number)
        self.new_Y = linear_interp(self.new_X)
        return self.new_X, self.new_Y


class chrom_plotter(diagram_base):
    def __init__(self, series, width=1800, height=800):
        super().__init__(width, height)
        self.selection_range = {}
        self.chrom_is_done = False
        self.column_volume = 7
        self.points_font_size = '14'
        self.selected_line = 'UV 1_260'
        self.inter_points_number = 240
        self.chrom_colors_list = ['lightgreen', 'orange', 'lightblue', 'yellow', 'red', 'white', 'pink', 'white']
        self.init_data = series
        self.interpolate_series()
        self.hiding_rect_after_up = True
        self.draw_axis()
        self.type_chrom_layer = self.image.add_layer()
        self.line_chrom_layer = self.image.add_layer()
        self.frac_chrom_layer = self.image.add_layer()
        self.axis_points_layer = self.image.add_layer()
        self.sel_line_layer = self.image.add_layer()
        self.set_type_chrom_colors()
        self.type_font_size = 20
        self.draw_type_chrom_layer()
        self.convert_series_data()
        if series != {}:
            self.draw_chrom_lines()
            self.draw_axis_points()

    def init_chrom(self, series):
        try:
            self.init_data = series
            self.interpolate_series()
            self.set_type_chrom_colors()
            self.draw_type_chrom_layer()
            self.convert_series_data()
            self.draw_chrom_lines()
            self.draw_axis_points()
            self.chrom_is_done = True
        except:
            self.chrom_is_done = False


    def interpolate_series(self):
        out_series = {}
        for key in self.init_data.keys():
            if key.find('Frac') == -1:
                k_keys = list(self.init_data[key])
                df = pd.DataFrame(
                    {
                        'x': self.init_data[key][k_keys[0]],
                        'y': self.init_data[key][k_keys[1]]
                    }
                )
                df['x'] = df['x'].astype(float)
                df['y'] = df['y'].astype(float)
                line_x, line_y = interpolate_crom_line(list(df['x']), list(df['y']),
                                                       point_numbers=self.inter_points_number)()
                out_series[key] = {
                    k_keys[0]: line_x,
                    k_keys[1]: line_y
                }
                self.init_data[key] = out_series[key].copy()



    def convert_series_data(self):
        self.lines = {}
        for key in self.init_data.keys():
            self.lines[key] = coord_scaler(self.init_data[key], self.draw_space)


    def set_type_chrom_colors(self):
        self.chrom_colors = {}
        for i, key in enumerate(self.init_data.keys()):
            self.chrom_colors[key] = self.chrom_colors_list[i]


    def draw_axis_points(self):
        k_keys = list(self.init_data[self.selected_line].keys())
        x = list(self.init_data[self.selected_line][k_keys[0]])
        y = list(self.init_data[self.selected_line][k_keys[1]])
        self.axis_points_layer.content = ''
        p_x = np.linspace(min(x), max(x), 30)
        p_y = np.linspace(min(y), max(y), 20)

        #scaler = coord_scaler(self.init_data, self.draw_space)
        coord_x, coord_y = set_draw_points(p_x, p_y, p_x, p_y, self.draw_space)

        self.x_axis_label = k_keys[0]
        self.y_axis_label = k_keys[1]
        self.draw_axis()

        for xx, c_x in zip(p_x, coord_x):
            self.type_chrom_layer.content += (f'<text x="{c_x}" '
                                              f'y="{self.draw_space["p2"][1] + self.axis_font_size}" '
                                              f'fill="{self.axis_line_color}" '
                                              f'font-size="{self.points_font_size}">{round(xx, 0)}</text>')

        for yy, c_y in zip(p_y, coord_y):
            self.type_chrom_layer.content += (f'<text x="{self.draw_space["p1"][0] - self.axis_font_size*2}" '
                                              f'y="{c_y + self.axis_font_size}" '
                                              f'fill="{self.axis_line_color}" '
                                              f'font-size="{self.points_font_size}">{round(yy, 0)}</text>')


    def draw_type_chrom_layer(self):
        self.type_chrom_layer.content = ''
        x = self.width // 3
        for key in self.init_data.keys():
            self.type_chrom_layer.content += (f'<text x="{x}" '
                                        f'y="{self.type_font_size + 10}" '
                                        f'fill="{self.chrom_colors[key]}" '
                                        f'font-size="{self.type_font_size}">{key}</text>')
            #x += len(key) * self.type_font_size
            x += 200

    def draw_fraction_layer(self, points):
        self.frac_chrom_layer.content = ''
        for x, y in points:
            self.frac_chrom_layer.content += (f'<line x1="{x}" '
                                            f'y1="{self.draw_space["p2"][1] - 2}" '
                                            f'x2="{x}" '
                                            f'y2="{self.draw_space["p2"][1] - self.height // 4}" '
                                            f'stroke="lightblue" '
                                            f'stroke-width="{self.axis_line_width}"/>')


    def draw_chrom_lines(self):
        self.line_chrom_layer.content = ''
        for key in self.lines.keys():
            if key.find('Frac') == -1:
                polyline = PolylineSVG(self.lines[key].points)
                polyline.color = self.chrom_colors[key]
                polyline.fill_color = 'none'
                if key == 'UV 1_260':
                    polyline.fill_color = self.chrom_colors[key]
                self.line_chrom_layer.content += polyline.svg_string()
            else:
                self.draw_fraction_layer(self.lines[key].points)


    def get_series_selected_points(self, key):
        series = self.init_data[key]
        k_keys = list(series.keys())
        points_ = {'x': [], 'y': []}
        points_['x'].append(self.selection_range['x1'])
        points_['y'].append(0.)
        for x, y in zip(series[k_keys[0]], series[k_keys[1]]):
            if x >= self.selection_range['x1'] and x <= self.selection_range['x2']:
                points_['x'].append(x)
                points_['y'].append(y)
        points_['x'].append(self.selection_range['x2'])
        points_['y'].append(0.)

        xx, yy = set_draw_points(np.array(points_['x']), np.array(points_['y']),
                                 np.array(points_['x']), series[k_keys[1]],
                                 self.selected_space)
        return [(i, j) for i, j in zip(xx, yy)]

    def draw_sel_line_layer(self):
        try:
            self.sel_line_layer.content = ''
            points = self.get_series_selected_points(self.selected_line)
            polyline = PolylineSVG(points)
            polyline.color = 'green'
            polyline.fill_color = 'yellow'
            self.sel_line_layer.content += polyline.svg_string()
            x, y = self.conc_b_selected['draw_point'][0], self.conc_b_selected['draw_point'][1]
            text = f'buff B: {self.conc_b_selected["conc_b"]}% '
            text += f'CV: {self.conc_b_selected["CV"]}'
            #print(x, y, text)
            self.sel_line_layer.content += (f'<text x="{x}" '
                                        f'y="{y}" '
                                        f'fill="{"white"}" '
                                        f'font-size="{self.type_font_size}">{text}</text>')
        except:
            pass


    def get_range_by_rect(self):
        x1, x2, y1, y2 = self.x_down_point, self.x_up_point, self.y_down_point, self.y_up_point
        series = self.init_data[self.selected_line]
        k_keys = list(series.keys())
        x_data = list(pd.Series(series[k_keys[0]]).astype(float))
        x_min, x_max = min(x_data), max(x_data)
        xcoord_min = self.draw_space['p1'][0]
        xcoord_max = self.draw_space['p2'][0]

        self.selected_space = {}
        self.selected_space['p1'] = (x1, self.draw_space['p1'][1])
        self.selected_space['p2'] = (x2, self.draw_space['p2'][1])

        x_cf = (x_max - x_min) / (xcoord_max - xcoord_min)
        xx_1 = x_cf * (x1 - xcoord_min) + x_min
        xx_2 = x_cf * (x2 - xcoord_min) + x_min
        self.selection_range = {}
        self.selection_range['x1'] = xx_1
        self.selection_range['x2'] = xx_2
        return xx_1, xx_2

    def culc_selected_conc_b(self):
        try:
            conc_b = self.init_data['Conc B']
            x1, x2 = self.selection_range['x1'], self.selection_range['x2']
            k_keys = list(conc_b.keys())
            conc_b = pd.DataFrame(conc_b)
            conc_b = conc_b[(conc_b[k_keys[0]] >= x1)&(conc_b[k_keys[0]] <= x2)]
            self.conc_b_selected = {}
            self.conc_b_selected['conc_b'] = float(round(conc_b[k_keys[1]].mean(), 1))
            self.conc_b_selected['CV'] = float(round(sum([x1, x2])/(self.column_volume * 2), 1))
            self.conc_b_selected['ml'] = float(round(sum([x1, x2])/2, 1))

            xx, yy = set_draw_points(np.array([self.conc_b_selected['ml']]), np.array([self.conc_b_selected['conc_b']]),
                                 np.array(self.init_data['Conc B'][k_keys[0]]),
                                 np.array(self.init_data['Conc B'][k_keys[1]]),
                                 self.draw_space)
            self.conc_b_selected['draw_point'] = (int(round(xx[0], 0)), int(round(yy[0], 0)))
        except:
            ui.notify('Необходимо добавить к хроматографии данные концентрации Б')


    def on_mouse_down(self, e):
        self.x_down_point = e.image_x
        self.y_down_point = e.image_y

    def on_mouse_move(self, e):
        x, y = e.image_x, e.image_y
        if self.down_btn:
            self.x_up_point = e.image_x
            self.y_up_point = e.image_y
            self.get_range_by_rect()
            self.draw_sel_line_layer()

    def on_mouse_up(self, e):
        self.x_up_point = e.image_x
        self.y_up_point = e.image_y
        self.get_range_by_rect()
        self.culc_selected_conc_b()
        self.draw_sel_line_layer()


    def get_data_to_base(self):
        export = {}
        export['chrom data'] = self.init_data.copy()
        export['select range'] = self.selection_range
        return export


class chrom_dialog():
    def __init__(self, rowdata):
        self.rowdata = rowdata
        self.data_from_base = {}

    #{'#': 19, "3'-end": 'BHQ2', "5'-end": 'Cy5', 'Amount, oe': '1-3', 'Lenght': '23', 'Name': 'IRF4-592-C',
    # 'Purification': 'Хроматография', 'Sequence': '[Alk]TAAAAGAAGGCAAATTCCCCTGT[BHQ2]', 'client id': 'НЦГИ',
    # 'input date': '09.10.2025', 'order id': 'НЦГИ/БЛМ/УТ1578', 'output date': '09.10.2025', 'status': 'in queue',
    # 'Order id': 7593, 'f_sequence': '[Cy5]TAAAAGAAGGCAAATTCCCCTGT[BHQ2]', 'CPG, mg': '5 mg',
    # 'Support type': 'bhq2_1000_hg', 'Do LCMS': True, 'Do synth': True, 'Do cart': False, 'Do hplc': True,
    # 'Do paag': False, 'Do sed': False, 'Do click': True, 'Do subl': True, 'Done LCMS': False, 'Done synth': True,
    # 'Done cart': False, 'Done hplc': False, 'Done paag': False, 'Done sed': False, 'Done click': False,
    # 'Done subl': False, 'Dens, oe/ml': 413, 'Vol, ml': '0.05', 'Purity, %': 50, 'Position': 'C3',
    # 'Purif type': 'Хроматография_Cy5', 'Date': '12.09.2025', 'Scale, OE': '1-3', 'Status': 'synthesis', 'DONE': False,
    # 'Wasted': False, 'Send': False, 'asm Sequence': '[10]TAAAAGAAGGCAAATTCCCCTGT', 'Synt number': '9', 'map #': 313.0,
    # 'Exist, oe': nan, 'sufficiency': nan, 'synt, positions': None}

    def show_content(self):
        print(self.data_from_base)
        with ui.dialog() as self.dialog:
            with ui.card().style('width: auto; max-width: none;'):
                #ui.textarea(value=f'{self.rowdata[0]}').style('width: 1000px')
                with ui.row():
                    ui.input(label='Название', value=self.rowdata[0]['Name']).style('width: 200px')
                    ui.input(label='ID', value=self.rowdata[0]['Order id']).style('width: 100px')
                    self.map_id = ui.input(label='MAP', value=self.rowdata[0]['map #']).style('width: 100px')
                    ui.input(label='Sequence', value=self.rowdata[0]['Sequence']).style('width: 400px')
                    ui.input(label='Lenght', value=self.rowdata[0]['Lenght']).style('width: 100px')
                    ui.input(label='Тип очистки', value=self.rowdata[0]['Purification']).style('width: 200px')
                    ui.input(label='Заказ', value=self.rowdata[0]['order id']).style('width: 200px')
                with ui.row():
                    self.file_input = ui.input(label='Название файла', value='').style('width: 500px')
                    self.column_input = ui.input(label='Объем колонки, мл', value='7',
                                                 on_change=self.on_column_volume_change,
                                                 ).style('width: 100px')
                with ui.row():
                    ui.upload(label='Загрузить хроматограмму',
                        on_upload=self.handle_upload).props("accept=.asc").classes("max-w-full")

                chrom_data = asc_reader('')
                self.chrom = chrom_plotter(chrom_data.get_series_by_contains(['UV', 'Conc']))

                with ui.row():
                    ui.button('Добавить', on_click=self.do_add_event)
                    ui.button('Отмена', on_click=self.dialog.close)

    def on_column_volume_change(self, e):
        try:
            if e.value != '':
                self.chrom.column_volume = int(e.value)
            else:
                self.chrom.column_volume = 7
        except:
            ui.notify('Нужно указать целочисленное значение объема колонки')

    def handle_upload(self, e: events.UploadEventArguments):
        self.filename = e.name
        self.file_input.value = self.filename
        text = e.content.readlines()
        data = []
        for row in text:
            data.append(str(row.decode()))
        chrom_data = asc_reader('')
        chrom_data.parse(data)
        self.chrom.init_chrom(chrom_data.get_series_by_contains(['UV', 'Conc']))

    def on_send_chrom_data(self, data):
        print(data)

    def do_add_event(self):
        data, out = {}, {}
        if self.chrom.chrom_is_done:
            if self.chrom.selection_range != {}:
                data['init_data'] = self.chrom.init_data.copy()

                series = {}
                for key in data['init_data'].keys():
                    series[key] = {}
                    for k in data['init_data'][key].keys():
                        series[key][k] = list(data['init_data'][key][k])

                data['init_data'] = series

                data['selection_range'] = self.chrom.selection_range
                data['conc_b_selected'] = self.chrom.conc_b_selected
                data['column_volume'] = self.chrom.column_volume
                data['filename'] = self.filename

                out['oligo_id'] = int(self.rowdata[0]['Order id'])
                if self.map_id.value == '':
                    ui.notify('MAP должен быть целым числом')
                    return False
                out['map_id'] = int(self.map_id.value)
                out['position'] = self.rowdata[0]['Position']
                out['chrom_data'] = data
            else:
                ui.notify('Нужно отметить диапазон элюции целевого продукта')
        self.on_send_chrom_data(out)
        self.dialog.close()
        return True


if __name__ in {"__main__", "__mp_main__"}:
    plot = diagram_base(1500, 500)
    plot.hiding_rect_after_up = True
    plot.draw_axis()

    #chrom_data = asc_reader('/home/alex/Documents/OligoSynt/Exports/RP_5ml_method_30_30_fractionation_1mll_alkine_new_050924 2_H5 050925 001.asc')
    chrom_data = asc_reader('/home/alex/Documents/OligoSynt/Exports/RP_5ml_method_30_30_fractionation_1mll_alkine_new_050924 10_B2_020925 001.asc')

    chrom = chrom_plotter(chrom_data.get_series_by_contains(['UV', 'Conc']))


    ui.run(port=8082)

