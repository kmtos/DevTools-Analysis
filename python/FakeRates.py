import os
import sys
import logging
import math

import ROOT

import operator

class FakeRates(object):
    '''Class to access the fakerates for a given lepton ID.'''

    def __init__(self,version):
        self.version = version
        self.fakehists = {'electrons':{},'muons':{},'taus':{}}
        self.fakehists_mc = {'electrons':{},'muons':{},'taus':{}}
        self.fakekey = '{num}_{denom}'
        # WZ fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        self.fake_rootfile = ROOT.TFile(fake_path)
        self.fakehists['electrons'][self.fakekey.format(num='WZMedium',denom='WZLoose')] = self.fake_rootfile.Get('e/medium/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='WZTight',denom='WZLoose')] = self.fake_rootfile.Get('e/tight/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='WZMedium',denom='WZLoose')] = self.fake_rootfile.Get('m/medium/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='WZTight',denom='WZLoose')] = self.fake_rootfile.Get('m/tight/fakeratePtEta')
        # H++ fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_hpp_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        if self.version=='80X':
            #fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_hpp_13TeV_Run2016B.root'.format(os.environ['CMSSW_BASE'])
            #fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_hpp_13TeV_Run2016BCD.root'.format(os.environ['CMSSW_BASE'])
            fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_dijet_hpp_13TeV_Run2016BCDEFGH.root'.format(os.environ['CMSSW_BASE'])
        self.fake_hpp_rootfile = ROOT.TFile(fake_path)
        self.fakehists['electrons'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_hpp_rootfile.Get('e/medium_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_hpp_rootfile.Get('e/tight_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_hpp_rootfile.Get('e/tight_medium/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_hpp_rootfile.Get('m/medium_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_hpp_rootfile.Get('m/tight_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_hpp_rootfile.Get('m/tight_medium/fakeratePtEta')
        # tau fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_w_tau_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        if self.version=='80X':
            fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_w_tau_13TeV_Run2016BCDEFGH.root'.format(os.environ['CMSSW_BASE'])
        self.fake_tau_rootfile = ROOT.TFile(fake_path)
        self.fakehists['taus'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLoose')] = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight',denom='HppLoose')] = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight',denom='HppMedium')] = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta_fromMC')

    def __exit__(self, type, value, traceback):
        self.__finish()

    def __del__(self):
        self.__finish()

    def __finish(self):
        self.fake_rootfile.Close()
        self.fake_hpp_rootfile.Close()
        self.fake_tau_rootfile.Close()

    def __get_fakerate(self,cand,num,denom):
        if cand.collName not in self.fakehists: return 0., 0.
        key = self.fakekey.format(num=num,denom=denom)
        if key not in self.fakehists[cand.collName]: return 0., 0.
        hist = self.fakehists[cand.collName][key]
        pt = cand.pt()
        eta = cand.eta()
        if pt > 100.: pt = 99.
        b = hist.FindBin(pt,abs(eta))
        return hist.GetBinContent(b), hist.GetBinError(b)

    def __get_fakerate_mc(self,cand,num,denom):
        if cand.collName not in self.fakehists_mc: return 0., 0.
        key = self.fakekey.format(num=num,denom=denom)
        if key not in self.fakehists_mc[cand.collName]: return 0., 0.
        hist = self.fakehists_mc[cand.collName][key]
        pt = cand.pt()
        eta = cand.eta()
        if pt > 100.: pt = 99.
        b = hist.FindBin(pt,abs(eta))
        return hist.GetBinContent(b), hist.GetBinError(b)

    def getFakeRate(self,cand,num,denom,doError=False):
        if not cand: return 0 # not defined
        val, err = self.__get_fakerate(cand,num,denom)
        return (val,min([1.,val+err]),max([0.,val-err])) if doError else val

    def getFakeRateMC(self,cand,num,denom,doError=False):
        if not cand: return 0 # not defined
        val, err = self.__get_fakerate_mc(cand,num,denom)
        return (val,min([1.,val+err]),max([0.,val-err])) if doError else val
