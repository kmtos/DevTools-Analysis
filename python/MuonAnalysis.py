# MuonAnalysis.py
import logging
import time
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from Candidates import *

import itertools
import operator

import ROOT

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
        self.addCandVar(label,'isTightMuon','isTightMuon','I')
        self.addCandVar(label,'isHighPtMuon','isHighPtMuon','I')
        self.addCandVar(label,'isPFMuon','isPFMuon','I')
        self.addCandVar(label,'isGlobalMuon','isGlobalMuon','I')
        self.addCandVar(label,'isTrackerMuon','isTrackerMuon','I')
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
        self.tree.add(lambda cands: cands[label].trackIso()/cands[label].pt(), '{0}_trackRelIso'.format(label), 'F')
