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


    def _getZ(self,genParticles):
        # first find leptons
        for lep1 in genParticles:
            if lep1.pdgId() not in [11,13,15]: continue
            for lep2 in genParticles:
                if lep2.pdgId() != -1 * lep1.pdgId(): continue

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
