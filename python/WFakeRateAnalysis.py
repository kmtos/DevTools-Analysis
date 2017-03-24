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

logger = logging.getLogger("WFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class WFakeRateAnalysis(AnalysisBase):
    '''
    Select a muon + lep for fake rate in a W control region
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','WFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','WFakeRateTree')
        super(WFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
            self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # lepton
        # mu tag
        self.addLeptonMet('wt')
        self.addLepton('t')
        self.tree.add(lambda cands: self.tightScale(cands['t']), 't_tightScale', 'F')

        # lep probe
        self.addLeptonMet('wl')
        self.addLepton('l')
        self.tree.add(lambda cands: self.passMedium(cands['l']), 'l_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['l']), 'l_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['l']), 'l_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['l']), 'l_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['l']), 'l_tightScale', 'F')

        # dilepton combination
        self.addDiLepton('z')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            't' : None,
            'l' : None,
            'z' : None,
            'wt': None,
            'wl': None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        muons = self.getCands(self.muons,self.passTight)
        elecs = self.getCands(self.electrons,self.passTight)
        if len(muons)<1 and len(elecs)<1: return candidate # need 1 tight muon/electron
        #if len(leps)!=2: return candidate # need 1 additional lep
        if len(leps)!=2: return candidate # need 1 additional lep

        # get invariant masses
        bestL = ()
        bestPt = 0
        for zpair in itertools.permutations(leps,2):
            if zpair[0].pt()<30: continue # trigger threshold
            if zpair[1].pt()<10: continue
            #if zpair[1].collName=='taus' and zpair[1].pt()<20: continue
            #if zpair[0].collName==zpair[1].collName: continue # me/em
            if zpair[0].charge()!=zpair[1].charge(): continue # same sign
            if not self.passTight(zpair[0]): continue # first must be tight
            z = DiCandidate(*zpair)
            if z.deltaR()<0.5: continue # well separated
            pt = zpair[1].pt()
            if pt>bestPt:
                bestL = zpair
                bestPt = pt

        if not bestL: return candidate # need a z candidate



        z = bestL
        candidate['t'] = z[0]
        candidate['l'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['wt'] = MetCompositeCandidate(self.pfmet,z[0])
        candidate['wl'] = MetCompositeCandidate(self.pfmet,z[1])


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
        #for coll in [self.muons,self.electrons,self.taus]:
        for coll in [self.muons,self.electrons]:
            cands += self.getCands(coll,passMode)
        return cands


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''.join([self.getCollectionString(cands[x]) for x in ['t','l']])
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
                'SingleElectron' : [
                    'Ele23_WPLoose_Gsf',
                ],
            }
        else:
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu24',
                    'IsoTkMu24',
                ],
                'SingleElectron' : [
                    'Ele27_WPTight_Gsf',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        if isinstance(cands['t'],Electron):
            datasets = ['SingleElectron']
        else:
            datasets = ['SingleMuon']
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands['t']]
        if isinstance(cands['t'],Electron):
            triggerList = ['Ele23_WPLoose'] if self.version=='76X' else ['Ele27Tight']
        else:
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
    parser.add_argument('--outputFile', type=str, default='wFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = WFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='WFakeRateTree',
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
