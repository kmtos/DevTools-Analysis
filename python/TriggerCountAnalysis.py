#!/usr/bin/env python
import argparse
import logging
import sys


from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight

from Candidates import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("TriggerCountAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class TriggerCountAnalysis(AnalysisBase):
    '''
    TriggerCount analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','triggerCountTree.root')
        outputTreeName = kwargs.pop('outputTreeName','TriggerCountTree')
        super(TriggerCountAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

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
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.DoubleMediumIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg', 'I')


    ############################
    ### select 3l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
        }

        return candidate

    ###########################
    ### analysis selections ###
    ###########################
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
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL',
                    'Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
                'SingleElectron' : [
                    'Ele23_WPLoose_Gsf',
                ],
                'Tau' : [
                    'DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg',
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
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
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
                'Tau' : [
                    'DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg',
                ],
            }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon', 
            'DoubleEG', 
            'MuonEG',
            'SingleMuon',
            'SingleElectron',
            'Tau',
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
        candList = [cands[c] for c in ['hpp1','hpp2','hm1']]
        numTaus = [c.collName for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['SingleMuSoup','SingleEleSoup','Mu17Mu8','Ele23Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','DoublePFTau35'] if self.version=='76X' else ['SingleMuSoup','SingleEleSoup','DoublePFTau35']
        else:
            triggerList = ['DoublePFTau35']
        return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)






def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('data'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='triggerCountTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    triggerCountAnalysis = TriggerCountAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='TriggerCountTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       triggerCountAnalysis.analyze()
       triggerCountAnalysis.finish()
    except KeyboardInterrupt:
       triggerCountAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
