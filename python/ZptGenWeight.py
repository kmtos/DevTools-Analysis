import sys
import os
import json

from math import floor

import ROOT


class ZptGenWeight(object):

    def __init__(self,version):
        self.version = version
        path = '{0}/src/DevTools/Analyzer/data/zpt_weights.root'.format(os.environ['CMSSW_BASE'])
        self.rootfile = ROOT.TFile(path)
        self.zpt_hist = self.rootfile.Get('zptmass_histo')

    def __exit__(self, type, value, traceback):
        self.__finish()

    def __del__(self):
        self.__finish()

    def __finish(self):
        self.rootfile.Close()

    def weight(self, genParticles):
        z = None
        w = 1.
        for gen in genParticles:
            if gen.pdgId()==23 and gen.numberOfDaughters()==2:
                z = gen
        if z:
            zpt = z.pt()
            zmass = z.mass()
            if zpt > 10000: zpt = 9999
            if zmass > 10000: zmass = 9999
            w = self.zpt_hist.GetBinContent(self.zpt_hist.FindBin(zmass,zpt))
        return w

    def weightNLO(self, genParticles):
        # from Qier
        z = None
        w = 1.
        for gen in genParticles:
            if gen.pdgId()==23 and gen.status()==62:
                z = gen

        if z:
            # note, should add gamma*
            zpt = z.pt()
            d = z.daughter_1()
            if abs(d)==11:
                if zpt<5:     w = 1.05
                elif zpt<8:   w = 1.10
                elif zpt<10:  w = 1.06
                elif zpt<15:  w = 1.03
                elif zpt<20:  w = 1.00
                elif zpt<30:  w = 0.92
                elif zpt<40:  w = 0.93
                elif zpt<50:  w = 0.94
                elif zpt<65:  w = 0.96
                elif zpt<85:  w = 1.00
                elif zpt<150: w = 1.02
            elif abs(d)==13:
                if zpt<5:     w = 1.06
                elif zpt<8:   w = 1.07
                elif zpt<10:  w = 1.05
                elif zpt<15:  w = 1.02
                elif zpt<20:  w = 0.97
                elif zpt<30:  w = 0.92
                elif zpt<40:  w = 0.90
                elif zpt<50:  w = 0.92
                elif zpt<65:  w = 0.94
                elif zpt<85:  w = 0.97
                elif zpt<150: w = 1.00

        return w
