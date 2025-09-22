import mzdatapy
import numpy as np
import json
from tqdm import tqdm
import pandas as pd
from lcms_chrom_data import interpolate_crom_line

def substract_bkg(data, bkg, treshold=3000, rt_min=50, rt_max=1500):
    ret = []

    mz_list = [i for i, f in enumerate(bkg) if f >= treshold]

    for d in tqdm(data):
        if not (int(round(d[1], 0)) in mz_list):
            if d[0] >= rt_min and d[0] <= rt_max:
                ret.append(d)

    return np.array(ret)

def get_intensity_map(data, low_treshold=1000, param=3): #mod
    max_mz = int(round(max(data[:, 1])*param, 0))
    max_t = int(round(max(data[:, 0]), 0))
    out = [[0. for t in range(0, max_t + 10)] for mz in range(0, max_mz + 10)]

    for d in data:
        if d[2] >= low_treshold:
            mz, t = int(round(d[1]*param, 0)), int(round(d[0], 0))
            out[mz][t] = d[2]
    return out

def find_neighbor(mz, t, map, param=3): #mod
    count = 0
    mz, t = int(round(mz*param, 0)), int(round(t, 0))
    for i in range(-1, 2):
        for j in range(-1, 2):
            if map[mz + i][t + j] > 0:
                count += 1
    return count * 100 / 9

def find_inner_points(data, map, neighbor_treshold=60, param=3):
    out = []

    for d in data:
        if find_neighbor(d[1], d[0], map, param=param) >= neighbor_treshold:
            out.append(d)

    return np.array(out)

class mzSpecDeconv():
    def __init__(self, mz_array, int_array, is_positive=False):
        self.data = pd.DataFrame({'mz': mz_array, 'intens': int_array,
                                  'class': np.zeros(mz_array.shape[0]),
                                  'mass': np.zeros(mz_array.shape[0])})
        self.is_positive = is_positive

    def __clusterize(self, data):

        clusters = []
        clusters.append([list(data.loc[0])])
        for index in range(1, data.shape[0]):
            finded = False
            for cl_index in range(len(clusters)):
                for mz_index in range(len(clusters[len(clusters) - cl_index - 1])):
                    mz_cl = clusters[len(clusters) - cl_index - 1][mz_index][0]
                    mz = data['mz'].loc[index]
                    if abs(mz - mz_cl) <= 1.05:
                        finded = True
                        clusters[len(clusters) - cl_index - 1].append(list(data.loc[index]))
                        data['class'].loc[index] = len(clusters) - cl_index - 1
                        break
                if finded:
                    break
            if not finded:
                clusters.append([list(data.loc[index])])
                data['class'].loc[index] = len(clusters) - 1

        return data

    def __clusterize_2(self, data):
        data['class'] = round(data['mz'], 0)
        return data

    def __compute_mass(self, data):

        classes = list(set(data['class']))
        data['charge'] = np.zeros(data['mz'].shape[0])
        data['mass'] = np.ones(data['mz'].shape[0])

        #print(data[data['class']==0.])

        if self.is_positive:
            sign = -1
        else:
            sign = 1

        for cl in classes:
            df = data[data['class'] == cl]
            if df.shape[0] > 3:
                df = df.sort_values(by='mz', ascending=False)
                #charge = round(1 / abs(df['mz'].values[0] - df['mz'].values[1]), 0)

                diff = pd.DataFrame(df['mz'])
                diff['diff'] = df['mz'].diff(periods=1)
                diff.dropna(inplace=True)
                diff = diff[diff['diff'] != 0]
                diff['charge'] = [abs(round(1/z, 0)) for z in diff['diff']]
                charge = diff['charge'].value_counts().idxmax()

                r_int = df['intens'] / df['intens'].sum()
                masses = df['mz'] * charge + sign * charge
                avg_mass = (masses * r_int).sum()

                data.loc[data['class'] == cl, 'charge'] = charge
                data.loc[data['class'] == cl, 'mass'] = avg_mass

        data['mono_mass'] = data['mz'] * data['charge'] + sign * data['charge']
        data = data[data['charge'] > 0]

        return data


    def deconvolute(self):
        self._data = self.data.sort_values(by='intens', ascending=False)
        self._data = self._data.reset_index()
        self._data = self._data.drop(['index'], axis=1)

        self._data = self.__clusterize(self._data)
        self._data = self.__compute_mass(self._data)

        #print(self._data)
        return self._data

    def deconvolute_2(self):
        self._data = self.data.sort_values(by='mz', ascending=True)
        self._data = self._data.reset_index()
        self._data = self._data.drop(['index'], axis=1)

        self._data = self.__clusterize_2(self._data)
        self._data = self.__compute_mass(self._data)

        #print(self._data)
        return self._data

    @staticmethod
    def drop_by_charge(data, max_charge=10):
        data = data[data['charge'] <= max_charge]
        return data

class oligosDeconvolution():
    def __init__(self, rt, mz, intens, is_positive=False, max_charge=10):
        rt = [round(i, 0) for i in rt]
        self.data = pd.DataFrame({'rt': rt, 'mz': mz, 'intens': intens})
        self.is_positive = is_positive
        self.max_charge = max_charge
        self.deconv_fast = True

    def deconvolute(self):
        rt_list = list(set(self.data['rt']))
        sum_data = np.array([])
        for rt in tqdm(rt_list, desc='Deconvolution:'):
            df = self.data[self.data['rt'] == rt]
            deconv = mzSpecDeconv(df['mz'], df['intens'], is_positive=self.is_positive)

            if self.deconv_fast:
                data = deconv.deconvolute_2()
            else:
                data = deconv.deconvolute()

            data = deconv.drop_by_charge(data, self.max_charge)
            data['rt'] = rt * np.ones(data['mz'].shape[0])

            if sum_data.shape[0] == 0:
                sum_data = data
            else:
                if data.shape[0] > 0:
                    sum_data = pd.concat([sum_data, data])

        sum_data = sum_data.sort_values(by='intens', ascending=False)
        sum_data = sum_data.reset_index()
        sum_data = sum_data.drop(['index'], axis=1)
        return sum_data

    @staticmethod
    def drop_data(data, mass_max, mass_min, rt_min, rt_max):
        df = data[data['mass'] <= mass_max]
        df = df[df['mass'] >= mass_min]
        df = df[df['rt'] <= rt_max]
        df = df[df['rt'] >= rt_min]
        return data.drop(list(df.index))

    @staticmethod
    def rt_filtration(data, rt_min, rt_max):
        df = data[data['rt'] <= rt_max]
        df = df[df['rt'] >= rt_min]
        return df

class zip_oligo_mzdata():
    def __init__(self, filename):
        self.filename = filename
        self.int_treshold = 1000
        self.max_mz = 3200
        self.rt_left = 0
        self.rt_right = 1500
        self.rt_interval = [self.rt_left, self.rt_right]
        self.bkg_treshold = 1000
        self.polish_param = 4
        self.delta_mass_treshold = 2
        self.neighbor_treshold = 60
        self.low_intens_treshold = 1000
        self.bkg_polish_count = 1

        self.rt_mul = 0.2
        self.mz_mul = 5
        self.mass_mul = 10
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
        self.compresed_data = self.compress_2(self.init_data)

    def from_base_file(self, filename):
        with open(filename, 'r') as file:
            data = file.read()
        self.base_data = json.loads(data)
        self.init_data = pd.DataFrame(self.base_data['init_data'])
        self.init_data = self.init_data.values


    def get_bkg(self, data):
        bkg = {}
        for i in range(int(round(self.max_mz, 0))):
            bkg[i] = 0
        for row in data:
            mz = int(round(row[1], 0))
            bkg[mz] += 1
        return list(bkg.values())


    def pipeline(self):
        out_dict = {}
        self.bkg = self.get_bkg(self.init_data)
        deconv_data = self.deconvolution()

        out_dict['init_zip'] = self.json_dumps_tuple_keys(self.compress_2(self.init_data))
        out_dict['deconv_zip'] = self.json_dumps_tuple_keys(self.compress_2(deconv_data[['rt', 'mass', 'intens']].values))
        out_dict['polish_zip'] = self.json_dumps_tuple_keys(self.compress_2(self.data))

        out_dict['init_chrom_line'] = json.dumps(self.get_chrom_data(self.init_data_to_df(self.init_data)))
        out_dict['polish_chrom_line'] = json.dumps(self.get_chrom_data(self.init_data_to_df(self.data)))
        out_dict['deconv_chrom_line'] = json.dumps(self.get_chrom_data(deconv_data))

        out_dict['oligo_ID'] = self.base_data['Order id']
        out_dict['map_ID'] = self.base_data['map id']
        out_dict['position'] = self.base_data['Position']

        return out_dict


    def json_dumps_tuple_keys(self, mapping):
        string_keys = {json.dumps(k): v for k, v in mapping.items()}
        return json.dumps(string_keys)

    def json_loads_tuple_keys(self, string):
        mapping = json.loads(string)
        return {tuple(json.loads(k)): v for k, v in mapping.items()}

    def compress_2(self, init_data):
        int_max = np.max(init_data[:, 2])

        rt = np.round(init_data[:, 0] * self.rt_mul, 0).astype(int)
        mz = np.round(init_data[:, 1] * self.mz_mul, 0).astype(int)
        intens = np.round(init_data[:, 2] * self.intens_mul / int_max, 0).astype(int)

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

    def fast_deconvolution(self, data):
        x_list, mz_list, intens_list, keys = [], [], [], []
        for key, value in zip(data.keys(), data.values()):
            x_list.append(key[0])
            x_list.append(key[1])
            mz_list.append(key[2] / self.mz_mul)
            mz_list.append(key[2] / self.mz_mul)
            intens_list.append(value * (key[1] - key[0]))
            intens_list.append(value * (key[1] - key[0]))
            keys.append(','.join([str(key[0]), str(key[1])]))
            keys.append(','.join([str(key[0]), str(key[1])]))
        deconv = oligosDeconvolution(x_list, mz_list, intens_list, keys, is_positive=False, max_charge=15)
        deconv.deconv_fast = True
        deconv_data = deconv.deconvolute()
        out = {}
        for key, mass, intens in zip(deconv_data['keys'], deconv_data['mass'], deconv_data['intens']):
            k = key.split(',')
            out[(int(k[0]), int(k[1]), int(round(mass * self.mz_mul, 0)))] = intens
        return out

    def polish(self, rt_interval,
                    bkg_treshold,
               neighbor_treshold,
               neighbor_repeats,
               int_treshold):

        self.rt_interval = rt_interval
        self.bkg_treshold = bkg_treshold
        self.neighbor_treshold = neighbor_treshold
        self.low_intens_treshold = int_treshold
        self.bkg_polish_count = neighbor_repeats

        self.substract_bkg()

        for i in range(self.bkg_polish_count):
            map = get_intensity_map(self.data, low_treshold=self.low_intens_treshold, param=self.polish_param)
            self.data = find_inner_points(self.data, map, neighbor_treshold=self.neighbor_treshold, param=self.polish_param)

    def substract_bkg(self):
        self.data = substract_bkg(self.init_data, self.bkg, treshold=self.bkg_treshold,
                                  rt_min=self.rt_interval[0], rt_max=self.rt_interval[1])

    def deconvolution(self):
        retention_interval = self.rt_interval
        bkg_treshold = self.bkg_treshold
        neighbor_treshold = self.neighbor_treshold
        bkg_polish_repeats = self.bkg_polish_count
        low_intens_treshold = self.low_intens_treshold
        self.polish(retention_interval, bkg_treshold, neighbor_treshold,
                              bkg_polish_repeats, low_intens_treshold)
        deconv = oligosDeconvolution(self.data[:, 0], self.data[:, 1], self.data[:, 2], is_positive=False, max_charge=10)
        deconv.deconv_fast = True
        deconv_data = deconv.deconvolute()
        return deconv_data

    def get_chrom_data(self, data, point_numbers=100):
        data = data.sort_values(by='rt', ascending=True)
        chrom_line = {'rt': [], 'tic': []}
        for rt in tqdm(data['rt'].unique()):
            df = data[data['rt'] == rt]
            chrom_line['rt'].append(rt)
            chrom_line['tic'].append(df['intens'].sum())
        min_value = min(chrom_line['tic'])
        chrom_line['tic'] = [i - min_value for i in chrom_line['tic']]
        model = interpolate_crom_line(chrom_line['rt'], chrom_line['tic'], point_numbers=point_numbers)()
        return {'rt': [round(float(i),2) for i in model[0]],
                'tic': [round(float(i),2) for i in model[1]]}

    def get_mz_zip_data(self, data):
        x_vals = []
        y_vals = []
        for key, value in zip(data.keys(), data.values()):
            x0, x1 = 0, value * (key[1] - key[0])
            y0, y1 = key[2] / self.mz_mul, key[2] / self.mz_mul
            x_vals.extend([x0, x1, None])  # None чтобы отделить линии
            y_vals.extend([y0, y1, None])
        return x_vals, y_vals

    def integrate_mz_data(self, data, a, b):
        target, total = 0., 0.
        for key, value in zip(data.keys(), data.values()):
            y = key[2] / self.mz_mul
            if (y >= a) and (y <= b):
                target += value
            total += value
        return target, total


    def init_data_to_df(self, data):
        return pd.DataFrame({
            'rt': data[:, 0],
            'mz': data[:, 1],
            'intens': data[:, 2]
        })


