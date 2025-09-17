from nicegui import app, ui
from invoce_page import api_db_interface
import pandas as pd
from OligoMap_utils import stock_db_data
from datetime import datetime
import requests
from io import BytesIO


class stock_data_page_model(api_db_interface):
    def __init__(self):
        IP = app.storage.general.get('db_IP')
        port = app.storage.general.get('db_port')
        super().__init__(IP, port)
        self.pincode = app.storage.user.get('pincode')

        self.init_page_content()

    def init_page_content(self):
        self.input_date()

        with ui.grid(columns=2).style("grid-template-columns: 1400px 1000px"):

            ui.label('Таблица отгруженной продукции:').style('font-size: 20px')
            ui.label('Таблица списаний сырья:').style('font-size: 20px')

            self.set_order_tab(1400, 1100)
            self.set_stock_tab(1300, 1100)
            with ui.row():
                ui.button('Показать продукцию', color='orange', on_click=self.on_show_products)
                ui.button('Скачать', color='green', on_click=self.on_save_products)

            with ui.row():
                ui.button('Показать списания', color='orange', on_click=self.on_show_write_off)
                ui.button('Скачать', color='green', on_click=self.on_save_rowmat)
                self.total_cost = ui.input(label='Итог:')


    def input_date(self):
        with ui.row():
            with ui.input('Начало перидоа').style('font-size: 20px') as self.start_date:
                with ui.menu().props('no-parent-event') as menu:
                    with ui.date().bind_value(self.start_date):
                        with ui.row().classes('justify-end'):
                            ui.button('Close', on_click=menu.close).props('flat')
                with self.start_date.add_slot('append'):
                    ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')

            with ui.input('Конец периода').style('font-size: 20px') as self.end_date:
                with ui.menu().props('no-parent-event') as menu:
                    with ui.date().bind_value(self.end_date):
                        with ui.row().classes('justify-end'):
                            ui.button('Close', on_click=menu.close).props('flat')
                with self.end_date.add_slot('append'):
                    ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')

    def get_finished_orders(self):
        url = f'{self.api_db_url}/get_orders_by_status/{self.db_name}/finished'
        ret = requests.get(url, headers=self.headers())
        if ret.status_code == 200:
            return ret.json()
        else:
            return []


    def on_show_products(self):
        data = self.get_finished_orders()

        start_date = datetime.strptime(self.start_date.value, "%Y-%m-%d")
        end_date = datetime.strptime(self.end_date.value, "%Y-%m-%d")

        df = pd.DataFrame(data)
        df['Date'] = pd.to_datetime(df['input date'], format='%m.%d.%Y')
        conditions = (df['Date'] >= start_date) & (df['Date'] <= end_date)
        columns = ["#", "Name", "5'-end", "Sequence", "3'-end", "Amount, oe", "Purification", "Lenght", "order id",
                   "client id", "input date", "output date"]
        df = df[conditions][columns]

        self.invoce_content_tab.options['rowData'] = df.to_dict('records')
        self.invoce_content_tab.update()


    def on_show_write_off(self):
        db = stock_db_data()
        df_out = pd.DataFrame(db.get_in_out_tab('output_tab'))
        df_tab = pd.DataFrame(db.get_total_tab_data())
        df_out['date'] = pd.to_datetime(df_out['Date'])

        start_date = datetime.strptime(self.start_date.value, '%Y-%m-%d')
        end_date = datetime.strptime(self.end_date.value, '%Y-%m-%d')
        df_out = df_out[(df_out['date'] >= start_date)&(df_out['date'] <= end_date)]
        df_out = df_out[['Name', 'Amount', 'Unicode', 'Date', 'Time', 'User', '#']]

        g_out_df = df_out.groupby('Unicode').agg({
            'Name': 'first',
            'Amount': 'sum',
            'User': 'first'
                                                  })
        g_out_df.reset_index(inplace=True)

        m_df = pd.merge(g_out_df, df_tab, how='inner', on='Unicode')
        df = m_df[['Name_y', 'Amount', 'units', 'producer', 'supplyer', 'price', 'Unicode']]
        df['sum'] = df['Amount'] * df['price']
        self.total_cost.value = str(df['sum'].sum())

        print(m_df.keys())

        self.ag_grid.options['rowData'] = df.to_dict('records')
        self.ag_grid.update()


    def set_order_tab(self, width=1200, height=1000):

        invoce_content_df = pd.DataFrame(
            {
                '#': [1],
                'Name': [''],
                "5'-end": [''],
                'Sequence': [''],
                "3'-end": [''],
                'Amount, oe': ['5-10'],
                'Exist, oe': [0.],
                'Purification': ['Cart'],
                'Lenght': [''],
                'status': ['in queue'],
                'input date': [''],
                'output date': [''],
                'client id': [''],
                'order id': [''],
                'sufficiency': [0.],
                'synt, positions': [''],
            }
        )

        columnDefs = [
            {"field": "#", 'filter': 'agTextColumnFilter', 'floatingFilter': True,
             "checkboxSelection": True, "headerCheckboxSelection": True,},
            {"field": "Name", "headerName": "Наименование",
             'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "5'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Sequence", "headerName": "Последовательность",
             'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "3'-end", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Amount, oe", "headerName": "Количество, ОЕ",
             'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "Exist, oe", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "Purification", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "Lenght", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "status", 'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "input date", "headerName": "Дата поступления",
             'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "output date", "headerName": "Дата готовности",
             'editable': True, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "order id", "headerName": "Номер заказа",
             'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "client id", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
           # {"field": "synt, positions", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
           # {"field": "sufficiency", 'filter': 'agTextColumnFilter', 'floatingFilter': True},
        ]

        self.invoce_content_tab = ui.aggrid(
            {
                'columnDefs': columnDefs,
                'rowData': invoce_content_df.to_dict('records'),
                'rowSelection': 'multiple',
                "pagination": True
            }
        ,
        theme='alpine-dark').style(f'width: {width}px; height: {height}px;')
        self.invoce_content_tab.auto_size_columns = True


    def set_stock_tab(self, width=120, height=1000):
        tab_df = pd.DataFrame({
            '#': [],
            'Name': [],
            "Amount": [],
            'Unicode': [],
            "Date": [],
            "Time": [],
            'User': []
        })

        colDefs = [
            #{"field": "#", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Name_y", "headerName": "Наименование",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Amount", "headerName": "Количество",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "Unicode",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "units", "headerName": "Ед. измерения",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "producer", "headerName": "Производитель",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "supplyer", "headerName": "Поставщик",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "price", "headerName": "Цена, за единицу",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            {"field": "sum", "headerName": "Сумма",
             'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "Date", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "Time", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
            #{"field": "User", 'editable': False, 'filter': 'agTextColumnFilter', 'floatingFilter': True},
        ]

        self.ag_grid = ui.aggrid(
                    {
                        'columnDefs': colDefs,
                        'rowData': tab_df.to_dict('records'),
                        'rowSelection': 'multiple',
                        "pagination": True,
                        # "enableRangeSelection": True,
                    }
                    ,
                    theme='alpine-dark').style(f'width: {width}px; height: {height}px')  # alpine  material  quartz  balham
        self.ag_grid.auto_size_columns = True

    def save_excel(self, filename, data):
        buffer = BytesIO()
        data.to_excel(buffer, index=False)
        buffer.seek(0)
        ui.download(buffer.read(), filename=f'{filename}.xlsx')

    def on_save_products(self):
        self.save_excel('products', pd.DataFrame(self.invoce_content_tab.options['rowData']))

    def on_save_rowmat(self):
        self.save_excel('rowmaterials', pd.DataFrame(self.ag_grid.options['rowData']))