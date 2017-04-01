#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import sys
import itertools
import operator

sys.argv.append('-b')
import ROOT
sys.argv.pop()

logger = logging.getLogger("WTauFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class WTauFakeRateAnalysis(AnalysisBase):
    '''
    Select a muon + tau for fake rate in a W control region
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','WTauFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','WTauFakeRateTree')
        super(WTauFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # lepton
        # mu tag
        self.addLeptonMet('wm')
        self.addLepton('m')
        self.tree.add(lambda cands: self.tightScale(cands['m']), 'm_tightScale', 'F')

        # tau probe
        self.addLeptonMet('wt')
        self.addLepton('t')
        self.tree.add(lambda cands: self.passMedium(cands['t']), 't_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['t']), 't_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['t']), 't_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['t']), 't_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['t']), 't_tightScale', 'F')

        # dilepton combination
        self.addDiLepton('z')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            'm' : None,
            't' : None,
            'z' : None,
            'wm': None,
            'wt': None,
            'met': self.pfmet,
        }

        # get leptons
        muons = self.getCands(self.muons,self.passTight)
        #loosemuons = self.getCands(self.muons,self.passLoose)
        taus = self.getCands(self.taus,self.passLoose)
        if len(muons)!=1: return candidate # need 1 tight muon
        #if len(loosemuons)>1: return candidate # dont allow another loose muon
        if len(taus)<1: return candidate # need at least 1 tau 


        # quality
        if muons[0].pt()<25: return candidate
        if taus[0].pt()<20: return candidate
        if muons[0].charge()!=taus[0].charge(): return candidate
        z = DiCandidate(muons[0],taus[0])
        if z.deltaR()<0.5: return candidate

        candidate['m'] = muons[0]
        candidate['t'] = taus[0]
        candidate['z'] = z
        candidate['wm'] = MetCompositeCandidate(self.pfmet,muons[0])
        candidate['wt'] = MetCompositeCandidate(self.pfmet,taus[0])


        return candidate


    ##################
    ### lepton IDs ###
    ##################
    def passLoose(self,cand):
        return passHppLoose(cand)

    def passMedium(self,cand):
        return passHppMedium(cand)

    def passTight(self,cand):
        return passHppTight(cand)

    def looseScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDLooseIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedVeto',cand)
        else:                           return 1.

    def mediumScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedMedium',cand)
        else:                           return 1.

    def tightScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedTight',cand)
        else:                           return 1.

    def getPassingCands(self,mode):
        if mode=='Loose':
            passMode = self.passLoose
        elif mode=='Medium':
            passMode = self.passMedium
        elif mode=='Tight':
            passMode = self.passTight
        else:
            return []
        cands = []
        for coll in [self.muons,self.taus]:
            cands += self.getCands(coll,passMode)
        return cands


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = 'mt'
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        isData = self.event.isData()>0.5
        if self.version=='76X':
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
            }
        else:
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu24',
                    'IsoTkMu24',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'SingleMuon',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands['m']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('SingleMuon'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='wTauFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = WTauFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='WTauFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       dyAnalysis.analyze()
       dyAnalysis.finish()
    except KeyboardInterrupt:
       dyAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

