from nicegui import ui
from raw_mat_widget import rawMatWidget
from raw_mat_widget import infoPanel

class rawmaterial_panel():
    def __init__(self, backend_model, ip_addr, port):
        self.backend_model = backend_model
        self.model = {}
        self.ip_addr, self.port = ip_addr, port
        self.pincode = ''


    def widgets_stack(self):
        wi_w = 200
        with ui.row().classes("w-full"):
            self.stack_group_solvents = ui.label('Solvents and reagents:').style('font-size: 30px;')

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000001', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000007', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000004', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000010', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000018', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000009', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000008', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000006', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000011', self.pincode)

        with ui.row().classes("w-full"):
            self.stack_group_solvents = ui.label('Phosphoramidites:').style('font-size: 30px;')

        wi_w = 200

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000120', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000121', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000122', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000123', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000129', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000133', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000130', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000124', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000125', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000126', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000127', self.pincode)

        with ui.row().classes("w-full"):
            self.stack_group_solvents = ui.label('Azides:').style('font-size: 30px;')

        wi_w = 200

        with ui.grid(columns=7).classes("w-full").style(f"grid-template-columns: "
                                                        f"{wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px {wi_w}px"):
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000166', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000167', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000168', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000170', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000171', self.pincode)
            rawMatWidget(self.ip_addr, self.port, 'INIT_BASE_CODE_OLIGO_LAB_0000172', self.pincode)

    def init_frontend(self):
        if self.pincode != '':
            with (((ui.grid(columns=2).classes("w-full").style(f"grid-template-columns: {1550}px {1200}px")))):
                with ui.column():
                    self.widgets_stack()
                with ui.column():
                    ui.label('Info panel:').style('font-size: 30px;')
                    self.info_panel = infoPanel(self.ip_addr, self.port, self.pincode)
                    with ui.row():
                        ui.label('Search by name:').style('font-size: 30px;')
                        self.search_text = ui.input(label='', placeholder='Search name',
                                                    on_change=self.info_panel.search_by_name).style('font-size: 30px;').classes('w-[700px]')
                        #self.on_search = ui.button('Search', color="#00a100").classes('w-[200px]')

                    #self.search_text.on_value_change = self.info_panel.search_by_name






