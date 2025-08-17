from nicegui import events, ui, app
from OligoMap_utils import api_db_interface
import requests
import pandas as pd

class infoPanel(api_db_interface):

    def __init__(self, api_IP, db_port, pincode):
        super().__init__(api_IP, db_port)
        self.pincode = pincode
        self._bkg = 'infopanel_bkg_2.png'

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.width = 1000
        self.height = 800
        self.search_top = 10

        self.search_color = 'orange'
        self.selection_color = 'yellow'
        self.border_color = 'orange'

        self.image = ui.interactive_image(f'/img/{self._bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])
        self.base_layer = self.image.add_layer()
        self.search_layer = self.image.add_layer()
        self.search_selection_layer = self.image.add_layer()
        self.draw()

        self.main_tab_data_df = pd.DataFrame(self.get_all_data_in_tab('total_tab'))


    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')


    def draw_search_layer(self, df):
        self.search_layer.content = ""
        dd = df.to_dict('records')
        x, y = 20, 10
        self.search_dict = {}
        for row in dd:
            self.search_dict[(x, y, self.width - 20, y + 35)] = row

            self.search_layer.content += (f'<text x="{20}" y="{y + 23}" '
                                          f'fill="white" font-size="20">{row[1]}</text>')

            self.search_layer.content += (f'<rect x={10} y={y} '
                                        f'width={self.width - 20} height={35} '
                                        f'fill="{self.search_color}" fill-opacity="0.2"'
                                        f'stroke="{self.search_color}" stroke-width="2"/>')

            y += 45

    def draw_search_selection(self, x, y):
        self.search_selection_layer.content = ""
        for key in self.search_dict.keys():
            if (y >= key[1]) and (y <= key[3]):
                self.selected_row = self.search_dict[key]

                self.search_selection_layer.content += (f'<rect x={10} y={key[1]} '
                                              f'width={self.width - 20} height={35} '
                                              f'fill="{self.selection_color}" fill-opacity="0.2"'
                                              f'stroke="{self.selection_color}" stroke-width="2"/>')


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'mousemove':
            if self.search_layer.content != "":
                self.draw_search_selection(x, y)

        if e.type == 'click':
            pass


    def search_by_name(self, e):
        self.search_dataframe = self.main_tab_data_df[self.main_tab_data_df[1].str.contains(e.value,
                                                                                            case=False, na=False)]
        if e.value != '':
            self.draw_search_layer(self.search_dataframe)
        else:
            self.search_layer.content = ""


    def get_all_data_in_tab(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/{tab_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()




class rawMatWidget(api_db_interface):

    def __init__(self, api_IP, db_port, unicode, pincode):
        super().__init__(api_IP, db_port)

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.unicode = unicode
        self.lable = 'raw material'

        self.pincode = pincode
        self.widget_bkg = 'widget_bkg_2.png'

        self.width = 150
        self.height = 100
        self.title_height = 30

        self.title_color = 'orange'
        self.border_color = 'orange'

        self.image = ui.interactive_image(f'/img/{self.widget_bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])
        self.base_layer = self.image.add_layer()

        self.info_layer = self.image.add_layer()

        self.get_info_from_base()
        self.draw()
        self.draw_info()


    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=15 ry=15 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')
        self.base_layer.content += (f'<rect x=6 y=6 width={self.width - 12} height={self.title_height} rx=10 ry=10 '
                                    f'fill="none" '
                                   f'stroke="{self.title_color}" stroke-width="2"/>')

    def draw_info(self):
        self.info_layer.content = ""

        name = self.info_data[0][1]

        self.info_layer.content += (f'<text x="{15}" y="{25}" '
                                                  f'fill="white" font-size="12">{name}</text>')

        fill_color = 'green'
        fill_width_div = 1
        if self.remain['exist'] <= self.info_data[0][5]:
            fill_color = 'red'
            fill_width_div = 2

        self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 25}" '
                                    f'fill="white" font-size="10">{self.info_data[0][3]}</text>')
        amount = f"{round(self.remain['exist'], 2)} / {self.info_data[0][5]}"
        self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 40}" '
                                    f'fill="white" font-size="10">{amount}</text>')

        self.info_layer.content += (f'<rect x=8 y={self.title_height + 12} '
                                    f'width={self.width // fill_width_div - 18} height={self.title_height + 20} '
                                    f'fill="{fill_color}" fill-opacity="0.4"'
                                    f'stroke="{fill_color}" stroke-width="2"/>')


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'click':
            self.get_info_from_base()

    def get_info_from_base(self):
        self.info_data = self.get_unicode_data_in_tab()
        self.remain = self.get_remaining_stock(self.unicode)
        #print(self.info_data)
        #print(self.remain)

    def get_remaining_stock(self, unicode):
        url = f"{self.api_db_url}/get_remaining_stock/{self.db_name}/{unicode}"
        input_ret = requests.get(url, headers=self.headers())
        try:
            return input_ret.json()
        except:
            return {'exist': 0.}

    def get_unicode_data_in_tab(self):
        url = f'{self.api_db_url}/get_keys_data/{self.db_name}/total_tab/unicode/{self.unicode}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()