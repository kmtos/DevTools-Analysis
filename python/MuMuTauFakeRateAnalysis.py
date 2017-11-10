#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MuMuTauFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class MuMuTauFakeRateAnalysis(AnalysisBase):
    '''
    MuMuTauFakeRate analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','muMuTauFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MuMuTauFakeRateTree')
        super(MuMuTauFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[0], 'triggerEfficiencyMC', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[1], 'triggerEfficiencyMCUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[2], 'triggerEfficiencyMCDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[0], 'triggerEfficiencyData', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[1], 'triggerEfficiencyDataUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[2], 'triggerEfficiencyDataDown', 'F')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('z1')
        self.addDetailedMuon('z1')
        self.addLepton('z2')
        self.addDetailedMuon('z2')

        # tau
        self.addDiLepton('mt')
        self.addLepton('t')
        self.addDetailedTau('t')
        self.addLepton('m')
        self.addDetailedMuon('m')

        # met
        self.addMet('met')


    ############################
    ### Additional functions ###
    ############################
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

    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        self.addCandVar(label,'isLooseMuon','isLooseMuon','I')
        self.addCandVar(label,'isMediumMuon','isMediumMuon','I')
        self.addCandVar(label,'isMediumMuonICHEP','isMediumMuonICHEP','I')
        self.addCandVar(label,'isTightMuon','isTightMuon','I')
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
        self.addCandVar(label,'relPFIsoDeltaBetaR04','relPFIsoDeltaBetaR04','F')
        self.tree.add(lambda cands: cands[label].trackIso()/cands[label].pt() if cands[label].pt() else 0., '{0}_trackRelIso'.format(label), 'F')

    def passMuon(self,cand):
        if cand.pt()<3: return False
        if not cand.isPFMuon(): return False
        if not (cand.isGlobalMuon() or cand.isTrackerMuon()): return False
        return True

    def passTau(self,cand):
        if cand.pt()<10: return False
        if not cand.decayModeFinding(): return False
        return True

    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            't' : None,
            'z' : None,
            'm': None,
            'mt': None,
            'met': self.pfmet,
        }

        # get leptons
        muons = [m for m in self.muons if self.passMuon(m)]
        taus = [t for t in self.taus if self.passTau(t)]
        leps = muons+taus
        if len(muons)<2: return candidate
        if len(taus)<1: return candidate


        # get the candidates
        zCand = []
        for mm in itertools.permutations(muons,2):
            m1 = mm[0]
            m2 = mm[1]
            # charge OS
            if m1.charge()==m2.charge(): continue
            # require lead m pt>25
            if m1.pt()<20: continue
            if m2.pt()<10: continue
            # make composites
            z = DiCandidate(m1,m2)
            if not zCand: zCand = mm
            better = True
            if zCand[0].pt()>m1.pt():
                better = False
            elif zCand[1].pt()>m2.pt():
                better = False
            if better: zCand = mm
                
        if not zCand: return candidate

        m1 = zCand[0] if zCand[0].pt()>zCand[1].pt() else zCand[1]
        m2 = zCand[1] if zCand[0].pt()>zCand[1].pt() else zCand[0]

        # highest pt tau with DR>0.8 from selected muons
        goodTaus = [t for t in taus if deltaR(t.eta(),t.phi(),m1.eta(),m1.phi())>0.8 and deltaR(t.eta(),t.phi(),m2.eta(),m2.phi())>0.8]
        if len(goodTaus)==0: return candidate

        t = goodTaus[0]

        # find if there is a muon nearby
        otherMuons = [m for m in muons if m!=m1 and m!=m2 and deltaR(m.eta(),m.phi(),t.eta(),t.phi())<0.8]
        otherMuonsByDR = sorted(otherMuons, key=lambda m: deltaR(m.eta(),m.phi(),t.eta(),t.phi()))
        nearM = otherMuonsByDR[0] if otherMuonsByDR else Candidate(None)

        z = DiCandidate(m1,m2)
        if z.M()>120: return candidate
        if z.M()<60: return candidate

        candidate['z1'] = m1
        candidate['z2'] = m2
        candidate['t'] = taus[0]
        candidate['z'] = DiCandidate(m1,m2)
        candidate['m'] = nearM
        candidate['mt'] = DiCandidate(nearM,t) if len(otherMuonsByDR) else DiCandidate(Candidate(None),Candidate(None))

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
            }
        else:
            triggerNames = {
                'DoubleMuon'     : [
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['z1','z2']]
        triggerList = ['Mu17_Mu8'] if self.version=='76X' else ['Mu17Mu8']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)







def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy',version='80XMuMuTauTauZSkim'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='muMuTauFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    muMuTauFakeRateAnalysis = MuMuTauFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MuMuTauFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       muMuTauFakeRateAnalysis.analyze()
       muMuTauFakeRateAnalysis.finish()
    except KeyboardInterrupt:
       muMuTauFakeRateAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
