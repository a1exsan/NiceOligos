from nicegui import events, ui, app
from OligoMap_utils import api_db_interface
import requests
import pandas as pd
import json
import datetime


class writeOff_dialog(api_db_interface):
    def __init__(self, api_IP, db_port, pincode, data, write_off=True):
        super().__init__(api_IP, db_port)

        self.write_off = write_off
        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.pincode = pincode
        self.data = data

        if self.write_off:
            text = 'Списание'
            title = 'Списать материал:'
            self.tab_name = 'output_tab'
        else:
            text = 'Поступление'
            title = 'Добавить материал:'
            self.tab_name = 'input_tab'

        if self.data != {}:
            self.unicode = self.data[0][2]

            if len(data) == 1:
                with ui.dialog() as self.dialog:
                    with ui.card():
                        ui.label(title)
                        ui.label(f'{self.data[0][1]}')
                        self.amount = ui.input()
                        with ui.row():
                            ui.button(text, on_click=self.on_write_off_click)
                            ui.button('Отмена', on_click=self.dialog.close)
            elif len(data) > 1:
                with ui.dialog() as self.dialog:
                    with ui.card():
                        ui.label(title)
                        #ui.label(f'{self.data[0][1]}')

                        df = pd.DataFrame(self.data)
                        rowData = pd.DataFrame({
                            'name': df[1],
                            'amount': [0. for i in range(df.shape[0])],
                            'unicode': df[2]
                        })

                        colDefs = [
                            {"field": "name", 'editable': True},
                            {"field": "amount", 'editable': True},
                        ]

                        self.ag_grid = ui.aggrid(
                            {
                                'columnDefs': colDefs,
                                'rowData': rowData.to_dict('records'),
                                'rowSelection': 'multiple',
                                "pagination": True,
                                # "enableRangeSelection": True,
                            }
                            ,
                            theme='alpine-dark').classes('h-[800px]')  # alpine  material  quartz  balham
                        self.ag_grid.auto_size_columns = True
                        self.tab_rowdata = rowData.to_dict('records')
                        self.ag_grid.on("cellValueChanged", self.update_cell_data)

                        with ui.row():
                            ui.button(text, on_click=self.on_write_off_click_list)
                            ui.button('Отмена', on_click=self.dialog.close)

    def update_cell_data(self, e):
        self.tab_rowdata[e.args["rowIndex"]] = e.args["data"]
        #print(self.tab_rowdata)

    def on_write_off_click_list(self):
        self.substruct_from_stock(self.tab_name, self.tab_rowdata)
        #print(self.tab_rowdata)
        self.dialog.close()
        ui.run_javascript('location.reload();')

    def on_write_off_click(self):
        self.substruct_from_stock(self.tab_name,
                                  [{'name': self.data[0][1], 'amount': self.amount.value, 'unicode': self.data[0][2]}])
        self.dialog.close()
        ui.run_javascript('location.reload();')

    def substruct_from_stock(self, tab_name, rowdata):
        for row in rowdata:
            if float(row['amount']) > 0:
                url = f"{self.api_db_url}/insert_data/{self.db_name}/{tab_name}"
                r = requests.post(url,
                json=json.dumps(
                        [
                        row['name'], row['unicode'], row['amount'],
                        datetime.datetime.now().date().strftime(self.strftime_format),
                        datetime.datetime.now().time().strftime(self.time_format),
                        #user_id
                        self.get_user_id()
                        ]
                    )
                , headers=self.headers())

    def get_user_id(self):
        url = f"{self.api_db_url}/get_keys_data/{self.db_name}/users/pin/{self.pincode}"
        r = requests.get(url, headers=self.headers())
        return r.json()[0][2]


class raw_mat_base_widget(api_db_interface):
    def __init__(self, api_IP, db_port, unicode, pincode):
        super().__init__(api_IP, db_port)

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.unicode = unicode
        self.pincode = pincode

    #[874, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-04', '12:02:17', '1783121115']
    #[879, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-10', '10:35:14', '1783121115']
    #[880, 'бензол хч', 'INIT_BASE_CODE_OLIGO_LAB_0000007', 1.0, '2025-03-10', '10:36:06', '1783121115']

    def get_info_from_base(self):
        self.info_data = self.get_unicode_data_in_tab()
        self.remain = self.get_remaining_stock(self.unicode)
        self.output_data = self.get_unicode_output_data(self.unicode)

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

    def get_unicode_output_data(self, unicode):
        #url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/output_tab'
        #ret = requests.get(url, headers=self.headers())
        #pd.DataFrame(ret.json()).to_csv('output_tab.csv', sep='\t')
        url = f'{self.api_db_url}/get_keys_data/{self.db_name}/output_tab/unicode/{unicode}'
        ret = requests.get(url, headers=self.headers())

        return ret.json()

    def group_remain_unicode_data(self):
        df = pd.DataFrame(self.output_data)
        if df.shape[0] > 0:
            df['x_day'] = pd.to_datetime(df[4], format=self.strftime_format)
            df_g = df.groupby(pd.Grouper(key='x_day', freq='ME')).sum()
            df_g.reset_index(inplace=True)
            return list(df_g['x_day']), list(df_g[3])
        else:
            return [], []


class event():
    def __init__(self):
        self.type = 'null'
        self.x = 0
        self.y = 0

class menuItem():
    def __init__(self, name, color, canvas, x=10, y=10, w=110, h=40):
        self.name = name
        self.color = color
        self.canvas = canvas.add_layer()

        self.pos_x = x
        self.pos_y = y
        self.width = w
        self.height = h

        self.draw(event())

    def check_coord_click(self, e):
        if (e.image_x >= self.pos_x) and (e.image_x <= self.width + self.pos_x):
            if (e.image_y >= self.pos_y) and (e.image_y <= self.height + self.pos_y):
                return True
            else:
                return False
        else:
            return False

    def draw(self, e):
        fillop = 0.2
        if e.type == 'mousedown':
            fillop = 0.6

        self.canvas.content = ""
        self.canvas.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + 25}" '
                                      f'fill="white" font-size="20">{self.name}</text>')
        self.canvas.content += (f'<rect x={self.pos_x} y={self.pos_y} '
                               f'width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="{self.color}" fill-opacity="{fillop}"'
                                   f'stroke="{self.color}" stroke-width="4"/>')

    def on_click(self, e):
        pass

    def do_mousedown(self, e):
        if self.check_coord_click(e):
            self.draw(e)
            self.on_click(e)

    def do_mouseup(self, e):
        self.draw(e)


class infoPanel_menu():
    def __init__(self, img):
        self.items = {}
        self.image = img
        self.items['show info'] = menuItem('Show info', 'teal', self.image)
        self.items['write-off'] = menuItem('Write-off', 'orange', self.image, x=150)
        self.items['write-in'] = menuItem('Write-in', 'teal', self.image, x=270)
        self.items['add-btn'] = menuItem('Add reagent', 'teal', self.image, x=410, w=130)
        self.items['edit-btn'] = menuItem('Edit reagent', 'orange', self.image, x=550, w=130)

    def do_mousedown(self, e):
        for item in self.items.values():
            item.do_mousedown(e)

    def do_mouseup(self, e):
        for item in self.items.values():
            item.do_mouseup(e)

# red (красный)
# orange (оранжевый)
# amber (янтарный)
# yellow (жёлтый)
# lime (лаймовый)
# green (зелёный)
# emerald (изумрудный)
# teal (сине-зелёный)
# cyan (голубой)
# sky (небесный)
# blue (синий)
# indigo (индиго)
# violet (фиолетовый)
# purple (пурпурный)
# fuchsia (фуксия)
# pink (розовый)
# rose (роза)
# slate (сланцевый)
# gray (серый)
# zinc (цинковый)
# neutral (нейтральный)
# stone (каменный)

#Каждый цвет имеет оттенки с номерами 50, 100, 200… до 950, где 50 — самый светлый, 950 — самый тёмный. Например, для синего это могут быть классы text-blue-500 или bg-blue-200.


class rawMatList(api_db_interface):
    def __init__(self, api_IP, db_port, pincode, image, label_obj):
        super().__init__(api_IP, db_port)
        self.image = image
        self.label_obj = label_obj
        self.pincode = pincode
        self.visible = True

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.width = 1000
        self.height = 800
        self.search_top = 10
        self.scrol_width = 20
        self.menu_height = 100

        self.search_color = 'orange'
        self.selection_color = 'yellow'
        self.border_color = 'orange'
        self.selected_list_color = 'green'

        self.search_layer = self.image.add_layer()
        self.search_selection_layer = self.image.add_layer()
        self.selected_list_layer = self.image.add_layer()

        self.scroll_layer = self.image.add_layer()
        self.scroll_len = 30
        self.scroll_pos = 2
        self.scroll_yy = 2
        self.scroll_down = False
        self.scroll_border = self.height - self.menu_height - self.scroll_len - 5

        self.pushed_obj = {}
        self.selected_list = {}
        self.selected_row = {}

        self.main_tab_data_df = pd.DataFrame(self.get_all_data_in_tab('total_tab'))
        self.search_dataframe = self.main_tab_data_df
        self.search_list_index = 0

        self.on_mouse_down = self.on_down_event
        self.on_mouse_click = self.on_click_event
        self.draw_scroll()

    def set_visible(self, value):
        self.visible = value
        if not value:
            d = {}
            d['self.search_layer'] = self.search_layer.content
            d['self.search_selection_layer'] = self.search_selection_layer.content
            d['self.selected_list_layer'] = self.selected_list_layer.content
            d['self.scroll_layer'] = self.scroll_layer.content

            self.search_layer.content = ""
            self.search_selection_layer.content = ""
            self.selected_list_layer.content = ""
            self.scroll_layer.content = ""

            app.storage.user['raw_list_dict'] = d
        else:
            if 'raw_list_dict' in list(app.storage.user.keys()):
                d = app.storage.user.get('raw_list_dict')
                self.search_layer.content = d['self.search_layer']
                self.search_selection_layer.content = d['self.search_selection_layer']
                self.selected_list_layer.content = d['self.selected_list_layer']
                self.scroll_layer.content = d['self.scroll_layer']



    def do_mousedown(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            if self.check_scroll_catching(x, y):
                self.scroll_down = True
                self.scroll_yy = y
            self.on_mouse_down(e)

    def do_mousemove(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            if self.scroll_down:
                self.move_scroll(y)
            if self.search_layer.content != "":
                self.draw_search_selection(x, y)


    def do_mouseup(self, e):
        if self.visible:
            x, y = e.image_x, e.image_y
            self.scroll_down = False

    def do_mouseclick(self, e):
        if self.visible:
            self.on_mouse_click(e)


    def check_seldown_event(self, key, y):
        for key in self.search_dict.keys():
            if (y >= key[1]) and (y <= key[3]):
                return True
        return False

    def check_selection_list_area(self, e):
        if e.image_y >= self.menu_height + 10:
            return True
        else:
            return False

    def on_down_event(self, e):
        if self.check_selection_list_area(e):
            if not e.ctrl:
                self.selected_list = {}
            if 'key' in list(self.selected_row.keys()):
                if self.check_seldown_event(self.selected_row['key'], e.image_y):
                    self.pushed_obj = self.selected_row.copy()
            if 1 in list(self.selected_row.keys()):
                self.label_obj.set_text(f"Info panel:  <<  {self.pushed_obj[1]}  >>")
                ui.notify(self.selected_row)

        self.draw_selected_list()

    def on_click_event(self, e):
        if 'key' in list(self.pushed_obj.keys()):
            self.selected_list[self.pushed_obj['key']] = self.pushed_obj.copy()
        if e.shift:
            pass
        self.draw_selected_list()

    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')
        self.draw_scroll()

    def draw_scroll(self):
        self.scroll_layer.content = ""

        self.scroll_layer.content = (f'<rect x={self.width - self.scrol_width} y={self.menu_height + self.scroll_pos} '
                                     f'width={15} height={self.scroll_len} '
                                     f' rx=10 ry=10 fill="none" '
                                     f'stroke="{self.border_color}" stroke-width="4"/>')

    def draw_search_layer(self, df):
        self.search_layer.content = ""
        dd = df.to_dict('records')
        x, y = 20, self.menu_height
        self.search_dict = {}
        for row in dd[self.search_list_index:]:
            d = row.copy()
            d['key'] = (x, y, self.width - 20, y + 35)
            self.search_dict[(x, y, self.width - 20, y + 35)] = d

            self.search_layer.content += (f'<text x="{20}" y="{y + 23}" '
                                          f'fill="white" font-size="20">{row[1]}</text>')

            self.search_layer.content += (f'<rect x={10} y={y} '
                                          f'width={self.width - 20 - self.scrol_width} height={35} '
                                          f'fill="{self.search_color}" fill-opacity="0.2"'
                                          f'stroke="{self.search_color}" stroke-width="2"/>')

            # print(y, self.height)
            y += 45
            if y > self.height:
                break

    def draw_selected_list(self):
        self.selected_list_layer.content = ""
        for row in self.selected_list.values():
            # print(row)
            key = row['key']
            self.selected_list_layer.content += (f'<rect x={10} y={key[1]} '
                                                 f'width={self.width - 20 - self.scrol_width} height={35} '
                                                 f'fill="{self.selected_list_color}" fill-opacity="0.2"'
                                                 f'stroke="{self.selected_list_color}" stroke-width="2"/>')

    def draw_search_selection(self, x, y):
        self.search_selection_layer.content = ""
        for key in self.search_dict.keys():
            if (y >= key[1]) and (y <= key[3]):
                self.selected_row = self.search_dict[key]
                self.search_selection_layer.content += (f'<rect x={10} y={key[1]} '
                                                        f'width={self.width - 20 - self.scrol_width} height={35} '
                                                        f'fill="{self.selection_color}" fill-opacity="0.2"'
                                                        f'stroke="{self.selection_color}" stroke-width="2"/>')

    def check_scroll_catching(self, x, y):
        if (x >= self.width - self.scrol_width) and (x <= self.width - self.scrol_width + 15):
            if (y >= self.menu_height + self.scroll_pos) and (
                    y <= self.menu_height + self.scroll_pos + self.scroll_len):
                return True
            else:
                return False
        else:
            return False

    def move_scroll(self, y):
        if (self.scroll_pos > 0) and (self.scroll_pos < self.scroll_border):
            self.scroll_pos = self.scroll_pos + y - self.scroll_yy
            if self.scroll_pos <= 0:
                self.scroll_pos = 2
            if self.scroll_pos >= self.scroll_border:
                self.scroll_pos = self.scroll_border - 1
            self.scroll_yy = y

        if self.search_dataframe.shape[0] >= 16:
            index = (self.search_dataframe.shape[0] - 16) * (self.scroll_pos - 2) / self.scroll_border
            index = int(round(index, 0))
        else:
            index = 0
        self.search_list_index = index
        self.draw_scroll()
        self.draw_search_layer(self.search_dataframe)

    def search_by_name(self, e):
        if e.value != '':
            self.search_dataframe = self.main_tab_data_df[self.main_tab_data_df[1].str.contains(e.value,
                                                                                                case=False, na=False)]
        else:
            self.search_dataframe = self.main_tab_data_df

        self.draw_search_layer(self.search_dataframe)


    def get_all_data_in_tab(self, tab_name):
        url = f'{self.api_db_url}/get_all_tab_data/{self.db_name}/{tab_name}'
        ret = requests.get(url, headers=self.headers())
        return ret.json()


class descriptionPanel(raw_mat_base_widget):
    def __init__(self, api_IP, db_port, unicode, pincode, canvas, pos_x=20, pos_y=110, width=960, height=660):
        super().__init__(api_IP, db_port, unicode, pincode)

        self.main_layer = canvas
        self.visible = True

        self.pos_x = pos_x
        self.pos_y = pos_y
        self.width= width
        self.height = height
        self.color = 'neutral'#'orange' # neutral
        self.chart_color = 'orange'
        self.bar_color = 'green'
        self.chart_height = 250

        self.get_info_from_base()
        self.x_data, self.y_remain = self.group_remain_unicode_data()
        #print(self.info_data)
        #print(self.remain)
        #[[1, 'Ацетонитрил. ХЧ', 'INIT_BASE_CODE_OLIGO_LAB_0000001', '1L flask',
        #  'Ацетонитрил. ХЧ; 80 ppm H2O; для синтеза (необходимо сушить над ситами 3 или 4 ангстрема) и ВЭЖХ', 20, 1]]
        #{'exist': 21.0}
        self.draw()

    def draw_remain_monthly_data(self):
        main_x, main_y, main_w, main_h = self.pos_x + 10, self.pos_y + 45, self.width - 20, self.chart_height

        number = 12

        max_y_list = max(self.y_remain[-number:])
        self.avп_cons = round(sum(self.y_remain[-number:]) / number, 2)
        self.summ_cons = round(sum(self.y_remain[-number:]), 2)
        max_y = 190

        y_list = [int(round(i * max_y / max_y_list, 0)) for i in self.y_remain[-number:]]

        coord_x, x_width = 40, self.width // (number * 2)
        for x, y, data in zip(self.x_data[-number:], y_list, self.y_remain[-number:]):
            month_text = x.strftime(format="%m")
            self.main_layer.content += (f'<rect x={coord_x} y={main_h - y + self.pos_y + 20} '
                                        f'width={x_width} height={y} '
                                        f' rx=10 ry=10 fill="{self.bar_color}" fill-opacity="{0.6}"'
                                        f'stroke="{self.bar_color}" stroke-width="4"/>')
            self.main_layer.content += (f'<text x="{coord_x}" y="{main_h + self.pos_y + 40}" '
                                        f'fill="white" font-size="20">{month_text} m</text>')
            self.main_layer.content += (f'<text x="{coord_x}" y="{main_h - y + self.pos_y + 10}" '
                                        f'fill="white" font-size="20">{round(data, 2)}</text>')
            coord_x += x_width * 2

        self.main_layer.content += (f'<text x="{40}" y="{main_y + 25}" '
                                    f'fill="white" font-size="20"> средний расход: {self.avп_cons}</text>')
        self.main_layer.content += (f'<text x="{40}" y="{main_y + 45}" '
                                    f'fill="white" font-size="20"> расход за год: {self.summ_cons}</text>')

        self.main_layer.content += (f'<rect x={main_x} y={main_y} '
                                    f'width={main_w} height={main_h} '
                                    f' rx=10 ry=10 fill="{self.chart_color}" fill-opacity="{0.1}"'
                                    f'stroke="{self.chart_color}" stroke-width="4"/>')

    def draw_desc_string(self, text_line, pos_y):
        self.main_layer.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + pos_y}" '
                                f'fill="white" font-size="20">{text_line}</text>')

    def draw_description(self):
        init_s = self.info_data[0][4]
        curs, str_curs = 0, 0
        text = ''
        pos_y = self.height * 2 // 3
        drawed = False
        while True:
            space = init_s.find(' ', curs)
            if space == -1:
                space = len(init_s)
            if len(init_s[str_curs:space]) <= 70:
                text += f' {init_s[curs: space]}'
                drawed = False
            else:
                text += f' {init_s[curs: space]}'
                self.draw_desc_string(text, pos_y)
                drawed = True
                str_curs = space + 1
                text = ''
                pos_y += 25
            curs = space + 1
            if curs >= len(init_s):
                if not drawed and text != '':
                    self.draw_desc_string(text, pos_y)
                break

    def draw(self):
        fillop = 0.2

        self.main_layer.content = ""
        self.draw_description()
        self.draw_remain_monthly_data()

        remain = self.remain['exist']
        units = self.info_data[0][3]
        self.main_layer.content += (f'<text x="{self.pos_x + 8}" y="{self.pos_y + 25}" '
                                f'fill="white" font-size="20">Остаток на складе: {round(remain, 0)} {units}</text>')
        self.main_layer.content += (f'<rect x={self.pos_x} y={self.pos_y} '
                                f'width={self.width} height={self.height} '
                                f' rx=10 ry=10 fill="{self.color}" fill-opacity="{fillop}"'
                                f'stroke="{self.color}" stroke-width="4"/>')

    def set_visible(self, value):
        self.visible = value
        if value:
            self.draw()
        else:
            self.main_layer.content = ""



class infoPanel(api_db_interface):

    def __init__(self, api_IP, db_port, pincode, lbl_obj):
        super().__init__(api_IP, db_port)
        self.api_IP, self.db_port = api_IP, db_port
        self.pincode = pincode
        self.label_obj = lbl_obj
        self._bkg = 'infopanel_bkg_2.png'

        self.db_name = 'stock_oligolab_5.db'
        self.strftime_format = "%Y-%m-%d"
        self.time_format = "%H:%M:%S"

        self.width = 1000
        self.height = 800
        self.search_top = 10
        self.scrol_width = 20
        self.menu_height = 100

        self.search_color = 'orange'
        self.selection_color = 'yellow'
        self.border_color = 'orange'
        self.selected_list_color = 'green'

        self.image = ui.interactive_image(f'/img/{self._bkg}', on_mouse=self.mouse_handler,  # cross='red',
                                          events=['mousedown', 'mousemove', 'mouseup', 'click'])

        self.base_layer = self.image.add_layer()
        self.info_layer = self.image.add_layer()

        self.draw()
        self.info_menu = infoPanel_menu(self.image)
        self.rawMat = rawMatList(api_IP, db_port, pincode, self.image, self.label_obj)

        self.info_menu.items['show info'].on_click = self.on_show_info_menu_click
        self.info_menu.items['write-off'].on_click = self.on_write_off_stock
        self.info_menu.items['write-in'].on_click = self.on_write_in_stock


    def draw(self):
        self.base_layer.content = (f'<rect x=0 y=0 width={self.width} height={self.height} '
                                   f' rx=10 ry=10 fill="none" '
                                   f'stroke="{self.border_color}" stroke-width="4"/>')

    def on_show_info_menu_click(self, e):
        self.rawMat.set_visible(not self.rawMat.visible)

        if self.rawMat.visible:
            self.rawMat.draw_scroll()

        if not self.rawMat.visible:
            self.info_menu.items['show info'].color = 'red'
            self.info_menu.items['show info'].name = 'close info'

            if 2 in list(self.rawMat.pushed_obj.keys()):
                unicode = self.rawMat.pushed_obj[2]
            else:
                unicode = ''
            self.descript_panel = descriptionPanel(self.api_IP, self.db_port, unicode, self.pincode, self.info_layer)
        else:
            if self.descript_panel != None:
                self.descript_panel.set_visible(False)
            self.info_menu.items['show info'].color = 'teal'
            self.info_menu.items['show info'].name = 'show info'

    def on_show_info_menu_widget_click(self, data):

        self.rawMat.pushed_obj = data[0]
        self.rawMat.label_obj.set_text(f"Info panel:  <<  {self.rawMat.pushed_obj[1]}  >>")
        ui.notify(self.rawMat.pushed_obj)

        self.rawMat.set_visible(False)

        if not self.rawMat.visible:
            self.info_menu.items['show info'].color = 'red'
            self.info_menu.items['show info'].name = 'close info'

            unicode = data[0][2]

            self.descript_panel = descriptionPanel(self.api_IP, self.db_port, unicode, self.pincode, self.info_layer)
        else:
            if self.descript_panel != None:
                self.descript_panel.set_visible(False)
            self.info_menu.items['show info'].color = 'teal'
            self.info_menu.items['show info'].name = 'show info'


    def on_write_off_stock(self, e):
        if len(self.rawMat.selected_list.keys()) > 1:
            data = [i for i in self.rawMat.selected_list.values()]
            woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, data)
            woff.dialog.open()
            self.rawMat.selected_list = {}
        else:
            if self.rawMat.pushed_obj != {}:
                woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, [self.rawMat.pushed_obj])
                woff.dialog.open()
                self.rawMat.pushed_obj = {}
            else:
                ui.notify('Выберете объект со склада')


    def on_write_in_stock(self, e):
        if len(self.rawMat.selected_list.keys()) > 1:
            data = [i for i in self.rawMat.selected_list.values()]
            woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, data, write_off=False)
            woff.dialog.open()
            self.rawMat.selected_list = {}
        else:
            if self.rawMat.pushed_obj != {}:
                woff = writeOff_dialog(self.api_IP, self.db_port, self.pincode, [self.rawMat.pushed_obj],
                                       write_off=False)
                woff.dialog.open()
                self.rawMat.pushed_obj = {}
            else:
                ui.notify('Выберете объект со склада')


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'mousedown':
            self.rawMat.do_mousedown(e)
            self.info_menu.do_mousedown(e)

        if e.type == 'mousemove':
            self.rawMat.do_mousemove(e)

        if e.type == 'mouseup':
            self.rawMat.do_mouseup(e)
            self.info_menu.do_mouseup(e)

        if e.type == 'click':
            self.rawMat.do_mouseclick(e)




class rawMatWidget(raw_mat_base_widget):

    def __init__(self, api_IP, db_port, unicode, pincode, on_click=lambda unicode: ui.notify(unicode), lbl='raw mat'):
        super().__init__(api_IP, db_port, unicode, pincode)

        self.label = lbl
        self.widget_bkg = 'widget_bkg_4.png'

        self.on_mouse_click = on_click

        self.width = 119
        self.height = 49
        self.title_height = 20

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
                                   f'stroke="{self.border_color}" stroke-width="2"/>')
        #self.base_layer.content += (f'<rect x=6 y=6 width={self.width - 12} height={self.title_height} rx=10 ry=10 '
        #                            f'fill="none" '
        #                           f'stroke="{self.title_color}" stroke-width="2"/>')

    def draw_info(self):
        self.info_layer.content = ""
        name = self.label
        self.info_layer.content += (f'<text x="{10}" y="{15}" '
                                                  f'fill="white" font-size="12">{name}</text>')
        fill_color = 'green'
        fill_const = 20
        fill_width_div = 1
        if self.remain['exist'] <= self.info_data[0][5]:
            fill_color = 'red'
            fill_width_div = 2
            fill_const = 0
        #self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 25}" '
        #                            f'fill="white" font-size="10">{self.info_data[0][3]}</text>')
        amount = f"{round(self.remain['exist'], 2)} / {self.info_data[0][5]}"
        self.info_layer.content += (f'<text x="{15}" y="{self.title_height + 15}" '
                                    f'fill="white" font-size="10">{amount}</text>')
        self.info_layer.content += (f'<rect x=8 y={self.title_height} '
                                    f'width={self.width // fill_width_div - fill_const} height={self.title_height} '
                                    f'fill="{fill_color}" fill-opacity="0.4"'
                                    f'stroke="{fill_color}" stroke-width="2"/>')


    def on_mouse_click(self, unicode):
        pass#print('click')


    def mouse_handler(self, e: events.MouseEventArguments):
        x, y = e.image_x, e.image_y

        if e.type == 'click':
            self.get_info_from_base()
            self.on_mouse_click(self.info_data)
