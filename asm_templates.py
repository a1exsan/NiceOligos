import json
import re

class asm2000_method_builder():
    def __init__(self):
        self.template_data = ''
        self.params = {}

    def open_file(self, fn):
        with open(fn) as file:
            self.template_data = file.read()
        if self.template_data != '':
            self.parse()

    def parse(self):
        pattern = r"\{\{(.*?)\}\}"
        matches = re.findall(pattern, self.template_data)
        mod_list = [match for match in matches]
        self.params = {}
        for mod in mod_list:
            value = '{' + mod + '}'
            key = '{{' + mod + '}}'
            self.params[key] = json.loads(value)

    def generate_method(self, fn):
        method = self.template_data
        for key in self.params.keys():
            d = self.params[key]
            val = d[list(d.keys())[0]]
            method = method.replace(key, val)
        with open(fn, 'w') as file:
            file.write(method)


def test():
    template = asm2000_method_builder()
    template.open_file('templates/autosampler_75ul_template.pr2')
    template.params['{{"debl_repeats":"5"}}'] = {'debl_repeats': '6'}
    template.params['{{"debl_dispense_repeat":"250"}}'] = {'debl_dispense_repeat': '300'}
    template.params['{{"wash_after_debl_dispense_2":"700"}}'] = {'wash_after_debl_dispense_2': '1000'}
    template.generate_method('templates/autosampler_75ul_method1.pr2')
    print(template.params)


if __name__ == '__main__':
    test()