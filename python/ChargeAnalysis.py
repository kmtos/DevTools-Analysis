#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from Candidates import *

import itertools
import operator

import ROOT


logger = logging.getLogger("ChargeAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class ChargeAnalysis(AnalysisBase):
    '''
    Charge analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','chargeTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ChargeTree')
        super(ChargeAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
            self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu22Pass(), 'pass_IsoMu22', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu22Pass(), 'pass_IsoTkMu22', 'I')
            self.tree.add(lambda cands: self.event.Mu45_eta2p1Pass(), 'pass_Mu45_eta2p1', 'I')
            self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.Ele25_eta2p1_WPTight_GsfPass(), 'pass_Ele25_eta2p1_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_eta2p1_WPLoose_GsfPass(), 'pass_Ele27_eta2p1_WPLoose_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele45_WPLoose_GsfPass(), 'pass_Ele45_WPLoose_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # z leptons
        self.addDiLepton('z')
        #self.addDiCandVar('z','z1','z2','mass_uncorrected','mass','F',uncorrected=True)
        self.addLepton('z1')
        self.tree.add(lambda cands: self.passMedium(cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z1']), 'z1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z1']), 'z1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z1']), 'z1_tightScale', 'F')
        self.addLepton('z2')
        self.tree.add(lambda cands: self.passMedium(cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z2']), 'z2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z2']), 'z2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z2']), 'z2_tightScale', 'F')

        # met
        self.addMet('met')

    ############################
    ### select DY candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            'z'  : None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        if len(leps)<2: return candidate # need at least 2 leptons

        # get invariant masses
        bestZ = ()
        bestst = 0
        for zpair in itertools.combinations(leps,2):
            if zpair[0].collName!=zpair[1].collName: continue # SF
            z = DiCandidate(*zpair)
            st = z.st()
            if st>bestst:
                bestZ = zpair
                bestst = st

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])

        return candidate

    #################
    ### lepton id ###
    #################
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
        for coll in [self.electrons,self.muons]:
            cands += self.getCands(coll,passMode)
        return cands


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def twoLoose(self,cands):
        return len(self.getPassingCands('Loose'))>=2

    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
        if self.version=='76X':
            triggerNames = {
                'DoubleMuon'     : [
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
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
                'DoubleMuon'     : [
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'SingleMuon'     : [
                    'IsoMu22',
                    'IsoTkMu22',
                    'Mu45_eta2p1',
                    'Mu50',
                ],
                'SingleElectron' : [
                    'Ele25_eta2p1_WPTight_Gsf',
                    'Ele27_WPTight_Gsf',
                    'Ele27_eta2p1_WPLoose_Gsf',
                    'Ele45_WPLoose_Gsf',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        if not cands['z1']: return False
        if isinstance(cands['z1'],Electron):
            datasets = [
                'DoubleEG',
                'SingleElectron',
            ]
        else:
            datasets = [
                'DoubleMuon',
                'SingleMuon',
            ]

        # reject triggers if they are in another dataset
        # looks for the dataset name in the filename
        # for MC it accepts all
        reject = True if self.event.isData()>0.5 else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0]: reject = False
            for trigger in triggerNames[dataset]:
                var = '{0}Pass'.format(trigger)
                passTrigger = getattr(self.event,var)()
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0]: break
        return False

    def triggerEfficiency(self,cands):
        candList = [cands[c] for c in ['z1','z2']]
        if isinstance(candList[0],Electron):
            triggerList = ['Ele23_WPLoose','Ele17_Ele12'] if self.version=='76X' else ['SingleEleSoup','Ele23Ele12']
        else:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Mu17_Mu8'] if self.version=='76X' else ['SingleMuSoup','Mu17Mu8']
        return self.triggerScales.getDataEfficiency(triggerList,candList)









def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='chargeTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = ChargeAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ChargeTree',
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

