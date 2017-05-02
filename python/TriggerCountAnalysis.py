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
        self.addTriggers()

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
        isData = self.event.isData()>0.5
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
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    'IsoMu24',
                    'IsoTkMu24',
                    #'Mu45_eta2p1',
                    #'Mu50',
                ],
                'SingleElectron' : [
                    'Ele27_WPTight_Gsf',
                ],
                'Tau' : [
                    'DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg',
                ],
            }
            if isData and self.event.run()>=281639:
                triggerNames['MuonEG'] = [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                    'Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ',
                ]

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
        return self.checkTrigger(*datasets,**triggerNames)


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
