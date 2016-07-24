#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
import logging
import time
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from Candidates import *

import itertools
import operator

import ROOT

logger = logging.getLogger("TauAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class TauAnalysis(AnalysisBase):
    '''
    all taus
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','tTree.root')
        outputTreeName = kwargs.pop('outputTreeName','TTree')
        super(TauAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup analysis tree

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('t')
        self.addDetailedTau('t')

        # met
        self.addMet('met')

    ####################################################
    ### override analyze to store after every lepton ###
    ####################################################
    def perRowAction(self):
        '''Per row action, can be overridden'''
        for tau in self.taus:
            cands = {
                't': tau,
                'w': MetCompositeCandidate(self.pfmet,tau),
                'met': self.pfmet,
                'event': self.event,
            }
            if tau.pt()<20: continue
            self.tree.fill(cands,allowDuplicates=True)

        self.eventsStored += 1

    #################
    ### detailed  ###
    #################
    def addDetailedTau(self,label):
        '''Add detailed variables'''
        self.addCandVar(label,'againstMuonLoose3','againstMuonLoose3','I')
        self.addCandVar(label,'againstMuonTight3','againstMuonTight3','I')
        self.addCandVar(label,'againstElectronVLooseMVA6','againstElectronVLooseMVA6','I')
        self.addCandVar(label,'againstElectronLooseMVA6','againstElectronLooseMVA6','I')
        self.addCandVar(label,'againstElectronMediumMVA6','againstElectronMediumMVA6','I')
        self.addCandVar(label,'againstElectronTightMVA6','againstElectronTightMVA6','I')
        self.addCandVar(label,'againstElectronVTightMVA6','againstElectronVTightMVA6','I')
        self.addCandVar(label,'decayModeFinding','decayModeFinding','I')
        self.addCandVar(label,'byIsolationMVArun2v1DBoldDMwLTraw','byIsolationMVArun2v1DBoldDMwLTraw','F')
        self.addCandVar(label,'byVLooseIsolationMVArun2v1DBoldDMwLT','byVLooseIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byLooseIsolationMVArun2v1DBoldDMwLT','byLooseIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byMediumIsolationMVArun2v1DBoldDMwLT','byMediumIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byTightIsolationMVArun2v1DBoldDMwLT','byTightIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byVTightIsolationMVArun2v1DBoldDMwLT','byVTightIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'decayModeFindingNewDMs','decayModeFindingNewDMs','I')
        self.addCandVar(label,'byIsolationMVArun2v1DBnewDMwLTraw','byIsolationMVArun2v1DBnewDMwLTraw','I')
        self.addCandVar(label,'byVLooseIsolationMVArun2v1DBnewDMwLT','byVLooseIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byLooseIsolationMVArun2v1DBnewDMwLT','byLooseIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byMediumIsolationMVArun2v1DBnewDMwLT','byMediumIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byTightIsolationMVArun2v1DBnewDMwLT','byTightIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byVTightIsolationMVArun2v1DBnewDMwLT','byVTightIsolationMVArun2v1DBnewDMwLT','I')

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='tTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    mAnalysis = TauAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='TTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       mAnalysis.analyze()
       mAnalysis.finish()
    except KeyboardInterrupt:
       mAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

