import mzdatapy
import numpy as np
import json
from tqdm import tqdm


class zip_oligo_mzdata():
    def __init__(self, filename):
        self.filename = filename
        self.int_treshold = 1000
        self.max_mz = 3200
        self.rt_left = 0
        self.rt_right = 1500
        self.rt_interval = [self.rt_left, self.rt_right]
        self.rt_mul = 0.2
        self.mz_mul = 5
        self.intens_mul = 100
        self.len_param = 20
        self.treshold = 1

    def open_file(self, filename):
        self.filename = filename
        spec = mzdatapy.mzdata(filename, from_string='')
        self.init_data, self.bkg = spec.mzdata2tab(int_treshold=self.int_treshold,
                                                   max_mz=self.max_mz,
                                                   rt_left=self.rt_interval[0])
        self.compresed_data = self.compress_2()

    def from_string(self, content):
        spec = mzdatapy.mzdata('', from_string=content)
        self.init_data, self.bkg = spec.mzdata2tab(int_treshold=self.int_treshold,
                                                   max_mz=self.max_mz,
                                                   rt_left=self.rt_interval[0])
        self.compresed_data = self.compress_2()

    def json_dumps_tuple_keys(self, mapping):
        string_keys = {json.dumps(k): v for k, v in mapping.items()}
        return json.dumps(string_keys)

    def json_loads_tuple_keys(self, string):
        mapping = json.loads(string)
        return {tuple(json.loads(k)): v for k, v in mapping.items()}

    def compress_2(self):
        int_max = np.max(self.init_data[:, 2])

        rt = np.round(self.init_data[:, 0] * self.rt_mul, 0).astype(int)
        mz = np.round(self.init_data[:, 1] * self.mz_mul, 0).astype(int)
        intens = np.round(self.init_data[:, 2] * self.intens_mul / int_max, 0).astype(int)

        d = {}
        for x, y, z in zip(rt, mz, intens):
            if z > 0:
                d[(int(x), int(y))] = int(z)
        d_keys = list(d.keys())

        new_d = {}
        for key in tqdm(d_keys):
            if d[key] > 0:
                intens_value = 0
                for i in range(self.len_param):
                    new_key = (key[0] + i, key[1])
                    if new_key in d_keys:
                        intens_value += d[new_key]
                        d[new_key] = 0
                    else:
                        if i > self.treshold:
                            new_d[(key[0], new_key[0], key[1])] = intens_value
                        break
        return new_d