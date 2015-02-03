#!/usr/bin/python

import sys, glob, re
import numpy
from numpy import array
from scipy.interpolate import UnivariateSpline
import os.path

line_pattern = re.compile(r'^ +\d')
field_separator = re.compile(' +')

# TODO: fix splines at endpoints

def read_cur_file(filename):
    infile = open(filename, 'r')
    heating = ([],[])
    cooling = ([],[])
    current = heating
    cool = False
    prev_temp = -300
    for line in infile:
        if line_pattern.match(line.rstrip()):
            temperature, mag_sus = \
                map(float, field_separator.split(line.lstrip())[0:2])
            if temperature < prev_temp - 0.5 and not cool:
                current = cooling
                cool = True
            current[0].append(temperature)
            current[1].append(mag_sus)
            prev_temp = temperature
    infile.close()
    cooling[0].reverse()
    cooling[1].reverse()
    heating = (heating[0][1:], heating[1][1:],)
    cooling = (cooling[0][1:], cooling[1][1:],)
    return (map(array,heating), map(array,cooling))
    
class Furnace:

    @staticmethod
    def extend_data((temps, mss)):
        tlist = temps.tolist()
        mlist = mss.tolist()
        tlist = [(tlist[0] - 20)] + [(tlist[0] - 10)] + \
            tlist + [tlist[-1] + 10] + [tlist[-1] + 20]
        mlist = [mlist[0]] + [mlist[0]] + \
            mlist + [mlist[-1]] + [mlist[-1]]
        return (array(tlist), array(mlist))

    def __init__(self, filename):
        heat, cool = read_cur_file(filename)
        (self.heat_data, self.cool_data) = \
            map(Furnace.extend_data, (heat, cool))
        self.heat_spline = UnivariateSpline(*self.heat_data, s = 5)
        self.cool_spline = UnivariateSpline(*self.cool_data, s = 5)

    @staticmethod
    def correct_with_spline((temps, mss), spline):
        mss_corrected = numpy.zeros_like(mss)
        for i in range(0, len(temps)):
            mss_corrected[i] = mss[i] - spline(temps[i])
        return (temps, mss_corrected)

    def correct(self, heating, cooling):
        return (Furnace.correct_with_spline(heating, self.heat_spline),
                Furnace.correct_with_spline(cooling, self.cool_spline))

class TdmsData:

    @staticmethod
    def shunt_up(values):
        'Move values up to ensure all >=0'
        if len(values)==0: return values
        minimum = min(values)
        if (minimum < 0): values = [v - minimum for v in values]
        return values

    @staticmethod
    def shunt((heat, cool), offset):
        heat_s = (heat[0], [m + offset for m in heat[1]])
        cool_s = (cool[0], [m + offset for m in cool[1]])
        return (heat_s, cool_s)

    def make_zero_at_700(self):
        'Correct values for a zero susceptibility at/near 700 degrees'
        print self.name, self.data.keys(), self.data[700][0][1][:5]
        offset = -min(self.data[700][0][1][-5:])
        new_data = {}
        for temp in self.data.keys():
            new_data[temp] = TdmsData.shunt(self.data[temp], offset)
        self.data = new_data

    @staticmethod
    def filename_to_temp(filename):
        'Convert a filename to a temperature'
        leafname = os.path.basename(filename)
        m = re.search(r'^(\d+)[AB]?\.CUR$', leafname)
        if m==None: return None
        return int(m.group(1))

    def correct_for_volume(self, data):
        scale = self.nom_vol / self.real_vol
        return [scale * datum for datum in data]

    def set_oom(self, new_oom):
        scale = 10. **(self.oom - new_oom)
        new_data = {}
        for (temp, (heating, cooling)) in self.data.items():
            heating2 = (heating[0], [ms*scale for ms in heating[1]])
            cooling2 = (cooling[0], [ms*scale for ms in cooling[1]])
            new_data[temp] = (heating2, cooling2)
        self.data = new_data
        self.oom = new_oom

    def read_files(self, sample_dir):
        cur_files = glob.glob(os.path.join(sample_dir, '*.CUR'))
        for filename in cur_files:
            temperature = TdmsData.filename_to_temp(filename)
            if temperature == None: continue
            (heating, cooling) = read_cur_file(filename)
            if self.furnace != None:
                heating, cooling = self.furnace.correct(heating, cooling)
            #heating = (heating[0], TdmsData.shunt_up(heating[1]))
            #cooling = (cooling[0], TdmsData.shunt_up(cooling[1]))
            heating = (heating[0], self.correct_for_volume(heating[1]))
            cooling = (cooling[0], self.correct_for_volume(cooling[1]))
            self.data[temperature] = (heating, cooling)
        # self.make_zero_at_700()
        self.samples.sort()

    def __init__(self, furnace, sample_dir, real_vol=0.25, nom_vol = 10.):
        self.oom = -6. # order of magnitude
        self.name = os.path.basename(sample_dir)
        self.furnace = furnace
        self.data = {}
        self.samples = []
        self.nom_vol = nom_vol
        self.real_vol = real_vol
        if (sample_dir != None): self.read_files(sample_dir)

    @staticmethod
    def chop_data((temps, mss), min_temp, max_temp):
        temps_out = []
        mss_out = []
        for i in range(0, len(temps)):
            temp = temps[i]
            if temp>=min_temp and temp<=max_temp:
                temps_out.append(temp)
                mss_out.append(mss[i])
        return (array(temps_out), array(mss_out))

    @staticmethod
    def linear_fit(xs, ys):
        fit = numpy.polyfit(xs, ys, 1)
        poly = numpy.poly1d(fit.tolist())
        model_ys = poly(xs)
        mean_y = numpy.mean(ys)
        sserr = numpy.sum((ys - model_ys)**2)
        sstot = numpy.sum((ys - mean_y)**2)
        rsquared = 1 - sserr / sstot
        return poly, rsquared

    def curie_paramag(self, cycle, min_temp, max_temp):
        ''' Estimates a Curie temperature by linear fit to 
        inverse susceptibility over a given range.'''
        (temps, mss) = TdmsData.chop_data(self.data[cycle][1], min_temp, max_temp)
        poly, rsquared = TdmsData.linear_fit(temps, 1./mss)
        curie = poly.r[0] # x axis intercept
        return (curie, rsquared, poly)

    def curie_inflection(self, cycle, min_temp, max_temp):
        ''' Estimates a Curie temperature by finding an inflection
        point on a Hopkinson peak.'''
        all_temps, all_mss = self.data[cycle][1]
        (temps, mss) = TdmsData.chop_data((all_temps, all_mss), min_temp, max_temp)
        spline = UnivariateSpline(all_temps, all_mss, s=.1)
        derivs = [ spline.derivatives(t)[2] for t in temps ]
        spline2 = UnivariateSpline(temps, derivs, s=3)
        return (spline2.roots()[0], spline)

    def alteration(self, cycle):
        return self.data[cycle][1][1][0] -  self.data[cycle][0][1][0]

    def alterations(self, cycles):
        return [self.alteration(cycle) for cycle in cycles]
