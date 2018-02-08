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
        self.fakehists['electrons'][self.fakekey.format(num='HppMedium',denom='HppLooseNew')] = self.fake_hpp_rootfile.Get('e/medium_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight', denom='HppLooseNew')] = self.fake_hpp_rootfile.Get('e/tight_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppMedium',denom='HppLoose')]    = self.fake_hpp_rootfile.Get('e/medium_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight', denom='HppLoose')]    = self.fake_hpp_rootfile.Get('e/tight_loose/fakeratePtEta')
        self.fakehists['electrons'][self.fakekey.format(num='HppTight', denom='HppMedium')]   = self.fake_hpp_rootfile.Get('e/tight_medium/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(    num='HppMedium',denom='HppLooseNew')] = self.fake_hpp_rootfile.Get('m/medium_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(    num='HppTight', denom='HppLooseNew')] = self.fake_hpp_rootfile.Get('m/tight_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(    num='HppMedium',denom='HppLoose')]    = self.fake_hpp_rootfile.Get('m/medium_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(    num='HppTight', denom='HppLoose')]    = self.fake_hpp_rootfile.Get('m/tight_loose/fakeratePtEta')
        self.fakehists['muons'][self.fakekey.format(    num='HppTight', denom='HppMedium')]   = self.fake_hpp_rootfile.Get('m/tight_medium/fakeratePtEta')
        # tau fakerates
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_w_tau_13TeV_Run2015D.root'.format(os.environ['CMSSW_BASE'])
        if self.version=='80X':
            # NOTE: W tau performs poorly, z tau much better
            #fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_w_tau_13TeV_Run2016BCDEFGH.root'.format(os.environ['CMSSW_BASE'])
            fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_z_tau_13TeV_Run2016BCDEFGH.root'.format(os.environ['CMSSW_BASE'])
        self.fake_tau_rootfile = ROOT.TFile(fake_path)
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseNew')]    = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseNew')]    = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEta')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseNew')]    = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseNew')]    = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEta_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseNewDM0')] = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM0')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseNewDM0')] = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM0')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseNewDM0')] = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM0_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseNewDM0')] = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM0_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseNewDM1')] = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM1')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseNewDM1')] = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM1')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseNewDM1')] = self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM1_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseNewDM1')] = self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM1_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseNewDM10')]= self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM10')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseNewDM10')]= self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM10')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseNewDM10')]= self.fake_tau_rootfile.Get('medium_newloose/fakeratePtEtaDM10_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseNewDM10')]= self.fake_tau_rootfile.Get('tight_newloose/fakeratePtEtaDM10_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLoose')]       = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLoose')]       = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppMedium')]      = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLoose')]       = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLoose')]       = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEta_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppMedium')]      = self.fake_tau_rootfile.Get('tight_medium/fakeratePtEta_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseDM0')]    = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM0')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseDM0')]    = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM0')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseDM0')]    = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM0_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseDM0')]    = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM0_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseDM1')]    = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM1')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseDM1')]    = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM1')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseDM1')]    = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM1_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseDM1')]    = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM1_fromMC')
        self.fakehists['taus'][self.fakekey.format(   num='HppMedium',denom='HppLooseDM10')]   = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM10')
        self.fakehists['taus'][self.fakekey.format(   num='HppTight', denom='HppLooseDM10')]   = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM10')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppMedium',denom='HppLooseDM10')]   = self.fake_tau_rootfile.Get('medium_loose/fakeratePtEtaDM10_fromMC')
        self.fakehists_mc['taus'][self.fakekey.format(num='HppTight', denom='HppLooseDM10')]   = self.fake_tau_rootfile.Get('tight_loose/fakeratePtEtaDM10_fromMC')

        # mmtt
        fake_path = '{0}/src/DevTools/Analyzer/data/fakerates_mmtt_tau_13TeV_Run2016BCDEFGH.root'.format(os.environ['CMSSW_BASE'])
        self.fake_mmtt_tau_rootfile = ROOT.TFile(fake_path)
        self.fakehists['taus'][self.fakekey.format(   num='HaaTight',denom='HaaLoose')] = self.fake_mmtt_tau_rootfile.Get('nearMuonMedium_nearMuon/fakeratePtEta')
        self.fakehists_mc['taus'][self.fakekey.format(num='HaaTight',denom='HaaLoose')] = self.fake_mmtt_tau_rootfile.Get('nearMuonMedium_nearMuon/fakeratePtEta_fromMC')

    def __exit__(self, type, value, traceback):
        self.__finish()

    def __del__(self):
        self.__finish()

    def __finish(self):
        self.fake_rootfile.Close()
        self.fake_hpp_rootfile.Close()
        self.fake_tau_rootfile.Close()
        self.fake_mmtt_tau_rootfile.Close()

    def __get_fakerate(self,cand,num,denom,doDM=False,hists=None):
        if not hists: hists = self.fakehists
        if cand.collName not in hists:
            logging.warning('{0} not in {1}'.format(cand.collName, hists.keys().__repr__()))
            return 0., 0.
        if doDM and cand.collName=='taus':
            dm = cand.decayMode()
            if dm in [0,5]:
                denom = denom+'DM0'
            elif dm in [1,6]:
                denom = denom+'DM1'
            elif dm in [10]:
                denom = denom+'DM10'
        key = self.fakekey.format(num=num,denom=denom)
        if key not in hists[cand.collName]:
            logging.warning('{0} not in {1}'.format(key, hists[cands.collName].keys().__repr__()))
            return 0., 0.
        hist = hists[cand.collName][key]
        pt = cand.pt()
        eta = cand.eta()
        if pt > 100.: pt = 99.
        b = hist.FindBin(pt,abs(eta))
        if not b:
            logging.warning('Failed to find bin {0}, {1}'.format(pt,eta))
        val,err = hist.GetBinContent(b), hist.GetBinError(b)
        if not val:
            logging.warning('Invalid fakerate found {}+/-{} for {} {} {} {}'.format(val,err,cand.collName,pt,eta,key))
        return val,err

    def __get_fakerate_mc(self,cand,num,denom,doDM=False):
        return self.__get_fakerate(cand,num,denom,doDM=doDM,hists=self.fakehists_mc)

    def getFakeRate(self,cand,num,denom,doError=False,doDM=False):
        if not cand: return 0 # not defined
        val, err = self.__get_fakerate(cand,num,denom,doDM=doDM)
        return (val,min([1.,val+err]),max([0.,val-err])) if doError else val

    def getFakeRateMC(self,cand,num,denom,doError=False,doDM=False):
        if not cand: return 0 # not defined
        val, err = self.__get_fakerate_mc(cand,num,denom,doDM=doDM)
        return (val,min([1.,val+err]),max([0.,val-err])) if doError else val
