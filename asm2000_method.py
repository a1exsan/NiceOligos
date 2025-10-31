from nicegui import ui, app
from asm_templates import asm2000_method_builder
from synthesis_method import synth_base

class asm2000_method_page_model():
    def __init__(self):
        self.method_base = synth_base()
        self.method_rowdata = self.method_base.get_all_rowdata()
        self.method = asm2000_method_builder()
        self.page_content()
        self.load_template()

    def page_content(self):
        self.field_dict = {}

        with ui.row():
            with ui.card():
                ui.label('Deblocking').style('font-size: 20px')
                with ui.row():
                    self.field_dict['debl_dispense'] = ui.input(label='First dispense, msec').style(
                        'width: 150px; font-size: 20px')
                    self.field_dict['debl_repeats'] = ui.input(label='Repeats').style('width: 100px; font-size: 20px')
                    self.field_dict['debl_dispense_repeat'] = ui.input(label='Dispense repeats, msec').style(
                        'width: 150px; font-size: 20px')
                with ui.row():
                    self.field_dict['wash_after_debl_1'] = ui.input(label='Wash_1 solution').style(
                        'width: 100px; font-size: 20px')
                    self.field_dict['wash_after_debl_dispense_1'] = ui.input(label='Wash_1 dispense, msec').style(
                        'width: 150px; font-size: 20px')
                    self.field_dict['wash_after_debl_2'] = ui.input(label='Wash_2 solution').style(
                        'width: 100px; font-size: 20px')
                    self.field_dict['wash_after_debl_dispense_2'] = ui.input(label='Wash_2 dispense, msec').style(
                        'width: 150px; font-size: 20px')
                with ui.row():
                    self.field_dict['debl_flow_rate'] = ui.input(label='flow rate, ul/sec', value=100).style(
                        'width: 100px; font-size: 20px')


    def load_template(self):
        self.method_rowdata = self.method_base.get_all_rowdata()
        self.method.template_data = self.method_rowdata[0]['template']
        self.method.parse()
        for value in self.method.params.values():
            for key in value.keys():
                if key in self.field_dict:
                    self.field_dict[key].value = value[key]

