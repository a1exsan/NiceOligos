import pandas as pd
from nicegui import events, ui, app
import random
from datetime import datetime

class PolylineSVG:
    def __init__(self, points):
        self.points = points
        self.color = 'gray'
        self.fill_color = "gray"
        self.stroke_width = 2
        self.fill_opacity = 0.3

    def svg_string(self):
        # Construct points string for SVG
        points_str = ' '.join(f'{x},{y}' for x, y in self.points)
        return f'''
                <polyline points="{points_str}" 
                          fill="{self.fill_color}"
                          fill-opacity="{self.fill_opacity}"
                          stroke="{self.color}" 
                          stroke-width="{self.stroke_width}"/>
            </svg>
        '''


class history_stat_data_sorter():
    def __init__(self, data):
        self.data = data

        self.sorting_mode = 'total month'
        self.sorting_mode = 'total week'
        self.sorting_mode = 'total week'
        self.sorting_mode = 'last Y month'
        self.sorting_mode = 'last Y week'
        self.sorting_mode = 'total quarter'

    def get_sorting_data(self, status='finished'):
        y_data = [int(row['data'][status]) for row in self.data]
        x_data = [row['Date'] for row in self.data]

        df = pd.DataFrame({
            'x': x_data,
            'y': y_data
        })
        df['x_day'] = pd.to_datetime(df['x'], format='%d.%m.%Y')
        df = df.sort_values('x_day', ascending=True)

        ff = 'M'
        if self.sorting_mode.find('quarter') > -1:
            ff = 'Q'
        elif self.sorting_mode.find('month') > -1:
            ff = 'M'
        elif self.sorting_mode.find('week') > -1:
            ff = 'W'

        df_g = df.groupby(pd.Grouper(key='x_day', freq=ff)).max()
        df_g['y_max'] = df_g['y']
        df_min = df.groupby(pd.Grouper(key='x_day', freq=ff)).min()
        df_g['y_min'] = df_min['y']
        df_sum = df.groupby(pd.Grouper(key='x_day', freq=ff)).sum()
        df_g['y_sum'] = df_sum['y']

        df_g.reset_index(inplace=True)
        df_g['y_sub'] = df_g['y_max'] - df_g['y_min']

        if self.sorting_mode.find('Y') > -1:
            now_Y = datetime.now().strftime('%Y')
            df_g = df_g[df_g['x_day'].dt.strftime('%Y') == now_Y]

        return df_g


class invoceChart():
    def __init__(self):
        self.pincode = ''
        self.chart_bkg = 'chart_bkg_3.png'

        self.width = 2000
        self.height = 200
        self.height_shift = 30
        self.r_rad = 8

        self.sorting_mode = 'total month'

        self.image = ui.interactive_image(f'/img/{self.chart_bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])
        self.base_layer = self.image.add_layer()
        self.test_chart = self.image.add_layer()
        self.finished_bar_layer = self.image.add_layer()
        self.fin_oligos_chart = self.image.add_layer()
        self.radio_btns = self.image.add_layer()

        self.mode_radio_btns = {}
        self.init_mode_radio_btns()

        self.draw()
        self.draw_test_data()
        self.draw_radio_buttons()

    def draw_stat_data(self, history_stat_data):
        self.stat_data = history_stat_data

        y_data = [int(row['data']['finished']) for row in self.stat_data]
        x_data = [row['Date'] for row in self.stat_data]

        df = pd.DataFrame({
            'x': x_data,
            'y': y_data
        })
        df['x'] = pd.to_datetime(df['x'], format='%d.%m.%Y')
        df = df.sort_values('x', ascending=True)
        df = df.groupby('x').agg('max')
        df.reset_index(inplace=True)

        self.draw_bar_diagram([i for i in range(df.shape[0])], list(df['y']))

        data_sorter = history_stat_data_sorter(self.stat_data)
        data_sorter.sorting_mode = self.get_mode_radio_btns()
        data = data_sorter.get_sorting_data(status='finished')

        x_data = list(data['x'])
        y_data = list(data['y_sub'])

        self.draw_date_bar_chart(x_data, y_data)

    def get_date_text(self, date):
        if self.sorting_mode.find('Y') > -1:
            return date[:-5]
        elif self.sorting_mode.find('week') > -1:
            return date[:-5]
        else:
            return date

    def draw_date_bar_chart(self, date_list, y_list):

        self.finished_bar_layer.content = ""

        x_list = [i for i in range(len(date_list))]

        y_data = self.normalisation(y_list, self.height_shift, self.height - self.height_shift)
        x_data = self.normalisation(x_list, 100, self.width - 100)

        width = self.width // len(date_list) // 2
        if width > 70:
            width = 70

        for x, y, date, fin in zip(x_data, y_data, date_list, y_list):
            self.finished_bar_layer.content += (f'<rect x={x - width} '
                                               f'y={self.height - y} '
                                               f'width={width} '
                                               f'height={y - self.height_shift} '
                                               f'fill="green" '
                                               f'fill-opacity="0.7"'
                                               f'stroke="lightgreen" '
                                               f'stroke-width="2"/>')

            self.finished_bar_layer.content += (f'<text x="{x - width}" y="{self.height - self.height_shift + 20}" '
                                                  f'fill="white" font-size="14">{self.get_date_text(date)}</text>')

            delta = 25
            if self.height - y - self.height_shift <= 20:
                delta = 45
            self.finished_bar_layer.content += (f'<text x="{x - width // 2 - 10}" y="{self.height - y - self.height_shift + delta}" '
                                                f'fill="white" font-size="14">{fin}</text>')


    def draw_test_data(self):
        y_data = [200 - random.randint(50, 190) for i in range(200)]
        x_data = [i * 9 + 50 for i in range(200)]

        for i in range(1, 200, 1):
            x0 = x_data[i - 1]
            y0 = y_data[i - 1]
            x1 = x_data[i]
            y1 = y_data[i]
            line = f'M {x0} {y0} ' +  f'L {x1} {y1} '
            self.test_chart.content += f'<path d="{line}" stroke="gray" stroke-width="4" fill="none" />'


    def normalisation(self, z_data, z_min, z_max):
        z_range = max(z_data) - min(z_data)
        k = (z_max - z_min) / z_range
        return [int(round((z - min(z_data)) * k  + z_min)) for z in z_data]


    def draw_bar_diagram(self, x_data, y_data):

        y_top = y_data[len(y_data) - 1]

        y_data = self.normalisation(y_data, self.height_shift, self.height - self.height_shift)
        x_data = self.normalisation(x_data, 100, self.width - 100)

        self.test_chart.content = ''
        self.fin_oligos_chart.content = ''

        points = []
        for i in range(len(y_data)):
            x, y = x_data[i], y_data[i]
            y= self.height - y
            points.append((x, y))

        x, y = x_data[len(y_data) - 1], 0
        y = self.height - y - self.height_shift
        points.append((x, y))

        x, y = x_data[0], 0
        y = self.height - y - self.height_shift
        points.append((x, y))

        x, y = x_data[0], y_data[0]
        y = self.height - y - self.height_shift
        points.append((x, y))

        polyline = PolylineSVG(points)
        self.fin_oligos_chart.content += polyline.svg_string()

        self.draw_top_of_the_summit(x_data[len(y_data) - 1], y_data[len(y_data) - 1],
                                    text=str(y_top))


    def draw_top_of_the_summit(self, x, y, text=''):
        points = []

        x1, y1 = x, y

        x, y = x, y
        y = self.height - y - self.height_shift / 3
        points.append((x, y))

        #print(x, y)

        x, y = x, y + 10
        points.append((x, y))

        #print(x, y)

        x, y = x + 20, y - 5
        points.append((x, y))

        #print(x, y)

        x, y = x1, y1
        y = self.height - y - self.height_shift / 3
        points.append((x, y))

        #print(x, y)

        polyline = PolylineSVG(points)
        polyline.fill_color = 'red'
        polyline.fill_opacity = '0.8'
        polyline.color = 'red'

        self.test_chart.content += polyline.svg_string()
        self.test_chart.content += (f'<text x="{x + 22}" y="{y}" '
                                                  f'fill="white" font-size="16">{text}</text>')



    def draw(self):
        self.base_layer.content = f'<rect x=0 y=0 width=2000 height=200 fill="none" stroke="orange" stroke-width="4"/>'

    def init_mode_radio_btns(self):
        x, y = 1920, 50
        self.mode_radio_btns[(x, y)] = {'1':'month', '2':True}
        x, y = 1920, 80
        self.mode_radio_btns[(x, y)] = {'1':'week', '2':False}
        x, y = 1920, 110
        self.mode_radio_btns[(x, y)] = {'1':'Y month', '2':False}
        x, y = 1920, 140
        self.mode_radio_btns[(x, y)] = {'1':'Y week', '2':False}
        x, y = 1920, 170
        self.mode_radio_btns[(x, y)] = {'1':'quarter', '2':False}


    def get_mode_radio_btns(self):
        for key, value in zip(self.mode_radio_btns.keys(), self.mode_radio_btns.values()):
            if value['2']:
                return value['1']


    def draw_radio_buttons(self):
        self.radio_btns.content = ""
        color = 'orange'
        fell_color = 'red'
        text = 'total month'
        for key, value in zip(self.mode_radio_btns.keys(), self.mode_radio_btns.values()):
            text = value['1']
            if value['2']:
                fell_color = 'gray'
            else:
                fell_color = 'none'
            self.radio_btns.content += (f'<circle cx="{key[0]}" cy="{key[1]}" r="{self.r_rad}" fill="{fell_color}" '
                                    f'stroke="{color}" stroke-width="2" />')

            self.radio_btns.content += (f'<text x="{key[0] + 20}" y="{key[1]}" '
                                                  f'fill="white" font-size="14">{text}</text>')

    def push_radio_btn(self, x, y, r):
        count = 0
        for key in self.mode_radio_btns.keys():
            num = (key[0] - x) ** 2 + (key[1] - y) ** 2
            if num <= r ** 2:
                count = 1
        if count == 1:
            for key in self.mode_radio_btns.keys():
                num = (key[0] - x)**2 + (key[1] - y)**2
                if num <= r**2:
                    self.mode_radio_btns[key]['2'] = True
                else:
                    self.mode_radio_btns[key]['2'] = False
        self.draw_radio_buttons()
        self.sorting_mode = self.get_mode_radio_btns()
        self.draw_stat_data(self.stat_data)

    def mouse_handler(self, e: events.MouseEventArguments):

        x, y = e.image_x, e.image_y

        if e.type == 'click':
            self.push_radio_btn(x, y, self.r_rad)
