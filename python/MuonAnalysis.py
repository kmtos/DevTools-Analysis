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

logger = logging.getLogger("MuonAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class MuonAnalysis(AnalysisBase):
    '''
    all s
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','mTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MTree')
        super(MuonAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup analysis tree

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('m')
        self.addDetailedMuon('m')

        # met
        self.addMet('met')

    ####################################################
    ### override analyze to store after every lepton ###
    ####################################################
    def perRowAction(self):
        '''Per row action, can be overridden'''
        for muon in self.muons:
            cands = {
                'm': muon,
                'w': MetCompositeCandidate(self.pfmet,muon),
                'met': self.pfmet,
                'event': self.event
            }
            if muon.pt()<10: continue
            self.tree.fill(cands,allowDuplicates=True)

        self.eventsStored += 1

    #####################
    ### detailed muon ###
    #####################
    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        self.addCandVar(label,'isLooseMuon','isLooseMuon','I')
        self.addCandVar(label,'isMediumMuon','isMediumMuon','I')
        self.addCandVar(label,'isMediumMuonICHEP','isMediumMuonICHEP','I')
        self.addCandVar(label,'isTightMuon','isTightMuon','I')
        self.addCandVar(label,'isHighPtMuon','isHighPtMuon','I')
        self.addCandVar(label,'isPFMuon','isPFMuon','I')
        self.addCandVar(label,'isGlobalMuon','isGlobalMuon','I')
        self.addCandVar(label,'isTrackerMuon','isTrackerMuon','I')
        self.addCandVar(label,'susyMVA','susyMVA','F')
        self.addCandVar(label,'isSUSYMVAPreselection','isSUSYMVAPreselection','I')
        self.addCandVar(label,'miniIsolation','miniIsolation','F')
        self.addCandVar(label,'miniIsolationCharged','miniIsolationCharged','F')
        self.addCandVar(label,'miniIsolationNeutral','miniIsolationNeutral','F')
        self.addCandVar(label,'miniIsolationPhoton','miniIsolationPhoton','F')
        self.addCandVar(label,'miniIsolationPileup','miniIsolationPileup','F')
        self.addCandVar(label,'muonBestTrackType','muonBestTrackType','I')
        self.addCandVar(label,'segmentCompatibility','segmentCompatibility','F')
        self.addCandVar(label,'isGoodMuon','isGoodMuon','I')
        self.addCandVar(label,'highPurityTrack','highPurityTrack','I')
        self.addCandVar(label,'matchedStations','matchedStations','I')
        self.addCandVar(label,'validMuonHits','validMuonHits','I')
        self.addCandVar(label,'normalizedChi2','normalizedChi2','F')
        self.addCandVar(label,'validPixelHits','validPixelHits','I')
        self.addCandVar(label,'trackerLayers','trackerLayers','I')
        self.addCandVar(label,'pixelLayers','pixelLayers','I')
        self.addCandVar(label,'validTrackerFraction','validTrackerFraction','F')
        self.addCandVar(label,'bestTrackPtError','bestTrackPtError','F')
        self.addCandVar(label,'bestTrackPt','bestTrackPt','F')
        self.addCandVar(label,'trackerStandaloneMatch','trackerStandaloneMatch','F')
        self.addCandVar(label,'trackKink','trackKink','F')
        self.addCandVar(label,'relPFIsoDeltaBetaR03','relPFIsoDeltaBetaR03','F')
        self.addCandVar(label,'relPFIsoDeltaBetaR04','relPFIsoDeltaBetaR04','F')
        self.tree.add(lambda cands: cands[label].trackIso()/cands[label].pt(), '{0}_trackRelIso'.format(label), 'F')
        self.addCandVar(label,'jetPtRatio','jetPtRatio','F')
        self.addCandVar(label,'jetPtRel','jetPtRel','F')
        self.addCandVar(label,'jetNumberOfChargedDaughters','jetNumberOfChargedDaughters','F')
        self.addCandVar(label,'jetBtagCSV','jetBtagCSV','F')


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='mTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    mAnalysis = MuonAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MTree',
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
