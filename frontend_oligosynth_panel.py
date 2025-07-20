from nicegui import ui
import copy as cp

import xwell_plate_unit as XWells

class oligosynth_panel():
    def __init__(self, backend_model):
        self.backend_model = backend_model
        self.model = {}

        with ui.grid(columns=3):
            ui.label('Cell11')
            ui.label('Cell12')
            ui.label('Cell13')

            ui.label('Cell21')
            self.xwells_obj = XWells.XWell_plate(self, self.backend_model)
            ui.label('Cell23')

        with ui.column():
            pass


    def get_model(self):
        self.model['xwells_obj'] = self.xwells_obj.get_copy()
        #print(self.model['xwells_obj'])
        return self.model

    def set_model(self, model):
        self.xwells_obj.set_copy(model['xwells_obj'])
        #print(model['xwells_obj'])