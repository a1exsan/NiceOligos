from nicegui import app, ui, events
from OligoMap_utils import api_db_interface
from lcms_zip import zip_oligo_mzdata
import base64
import mzdatapy
import plotly.graph_objects as go


class lcms_analyser():
    def __init__(self):
        #IP = app.storage.general.get('db_IP')
        #port = app.storage.general.get('db_port')
        #super().__init__(IP, port)
        #self.pincode = app.storage.user.get('pincode')

        self.init_ui()

    def init_ui(self):
        with ui.row():
            self.file_input = ui.input(label='filename')
            ui.upload(label='Загрузить хроматограмму',
                  on_upload=self.handle_upload).props("accept=.mzdata.xml").classes("max-w-full")
            ui.button('show file', on_click=self.on_compress)
            ui.button('plot init data', on_click=self.on_plot_init_data)
        self.fig = go.Figure()
        self.fig.update_layout(
            template='plotly_dark',  # темная тема от Plotly
            plot_bgcolor='rgba(0,0,0,0)',  # прозрачный фон графика
            paper_bgcolor='black',  # фон всего графика
            font=dict(color='white')  # белый цвет текста и осей
        )
        self.fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
        self.plot = ui.plotly(self.fig).style('width: 1200px; height: 800px;')

    def on_plot_init_data(self):
        self.draw_lcms_zip_data(self.zip_data)

    def draw_lcms_zip_data(self, data):
        #for key, value in zip(data.keys(), data.values()):
        #    plt.plot([key[0] * 5, key[1] * 5], [key[2]/5, key[2]/5], '-', c='blue')
        for key, value in zip(data.keys(), data.values()):
            x0, x1 = key[0] / self.zip_lcms.rt_mul, key[1] / self.zip_lcms.rt_mul
            y0, y1 = key[2] / self.zip_lcms.mz_mul, key[2] / self.zip_lcms.mz_mul
            self.fig.add_trace(go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode='lines',
                line=dict(width=2, color='cyan  '),
                showlegend=False
            ))
        self.plot.update()



    def handle_upload(self, e: events.UploadEventArguments):
        self.file_input.value = e.name
        data = e.content.read()
        b64_string = base64.b64encode(data).decode()
        self.zip_lcms = zip_oligo_mzdata('')
        self.zip_lcms.from_string(b64_string)

    def on_compress(self):
        self.zip_data = self.zip_lcms.compress_2()
        print(self.zip_data)


if __name__ in {"__main__", "__mp_main__"}:
    ui.dark_mode(True)

    lcms = lcms_analyser()

    ui.run(port=8081)