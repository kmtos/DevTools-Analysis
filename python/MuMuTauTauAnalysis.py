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

logger = logging.getLogger("MuMuTauTauAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class MuMuTauTauAnalysis(AnalysisBase):
    '''
    MuMuTauTau analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','muMuTauTauTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MuMuTauTauTree')
        super(MuMuTauTauAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')

        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[0], 'triggerEfficiencyMC', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[1], 'triggerEfficiencyMCUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[2], 'triggerEfficiencyMCDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[0], 'triggerEfficiencyData', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[1], 'triggerEfficiencyDataUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[2], 'triggerEfficiencyDataDown', 'F')

        # h
        self.addComposite('h')
        self.addCompositeMet('hmet')
        self.tree.add(lambda cands: cands['hmet'].Mcat(2,3), 'hmet_mcat', 'F')

        # hpp leptons
        self.addDiLepton('amm')
        self.addCompositeMet('ammmet')
        self.addLepton('am1')
        self.addDetailedMuon('am1')
        self.addLepton('am2')
        self.addDetailedMuon('am2')

        # hmm leptons
        self.addDiLepton('att')
        self.addCompositeMet('attmet')
        self.tree.add(lambda cands: cands['attmet'].Mcat(0,1), 'attmet_mcat', 'F')
        self.addLepton('atm')
        self.addDetailedMuon('atm')
        self.addLepton('ath')
        self.addDetailedTau('ath')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

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
        self.tree.add(lambda cands: cands[label].trackIso()/cands[label].pt(), '{0}_trackRelIso'.format(label), 'F')

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
            'am1' : None,
            'am2' : None,
            'atm' : None,
            'ath' : None,
            'amm' : None,
            'ammmet' : None,
            'att' : None,
            'attmet' : None,
            'h' : None,
            'hmet' : None,
            'met': self.pfmet,
            'cleanJets' : [],
        }

        # get leptons
        muons = [m for m in self.muons if self.passMuon(m)]
        taus = [t for t in self.taus if self.passTau(t)]
        leps = muons+taus
        if len(muons)<3: return candidate
        if len(taus)<1: return candidate


        # get the candidates
        hCand = []
        mmDeltaR = 999
        ttDeltaR = 999
        for quad in itertools.permutations(leps,4):
            # require mmmt
            if not quad[0].__class__.__name__=='Muon': continue
            if not quad[1].__class__.__name__=='Muon': continue
            if not quad[2].__class__.__name__=='Muon': continue
            if not quad[3].__class__.__name__=='Tau': continue
            # charge OS
            if quad[0].charge()==quad[1].charge(): continue
            if quad[2].charge()==quad[3].charge(): continue
            # require lead m pt>25
            if quad[0].pt()<25: continue
            # make composites
            amm = DiCandidate(quad[0],quad[1])
            att = DiCandidate(quad[2],quad[3])
            # skim level selections
            #if amm.M()>30: continue
            #if amm.deltaR()>1.5: continue
            # choose best
            if not hCand: hCand = quad
            better = True
            if amm.deltaR()>mmDeltaR:
                better = False
            elif att.deltaR()>ttDeltaR:
                better = False
            if better: hCand = quad
                
        if not hCand: return candidate

        am1 = hCand[0] if hCand[0].pt()>hCand[1].pt() else hCand[1]
        am2 = hCand[1] if hCand[0].pt()>hCand[1].pt() else hCand[0]
        atm = hCand[2]
        ath = hCand[3]

        amm = DiCandidate(am1,am2)
        if amm.M()>30: return candidate
        #if amm.deltaR()>1.5: return candidate

        candidate['am1'] = am1
        candidate['am2'] = am2
        candidate['atm'] = atm
        candidate['ath'] = ath
        candidate['amm'] = DiCandidate(am1,am2)
        candidate['ammmet'] = MetCompositeCandidate(self.pfmet,am1,am2)
        candidate['att'] = DiCandidate(atm,ath)
        candidate['attmet'] = MetCompositeCandidate(self.pfmet,atm,ath)
        candidate['h'] = CompositeCandidate(am1,am2,atm,ath)
        candidate['hmet'] = MetCompositeCandidate(self.pfmet,am1,am2,atm,ath)

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,[am1,am2,atm,ath],0.4)

        return candidate


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
                    #'Mu45_eta2p1',
                    #'Mu50',
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
        candList = [cands[c] for c in ['am1','am2','atm','ath']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)







def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('haa',version='80XMuMuTauTau'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='muMuTauTauTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    muMuTauTauAnalysis = MuMuTauTauAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MuMuTauTauTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       muMuTauTauAnalysis.analyze()
       muMuTauTauAnalysis.finish()
    except KeyboardInterrupt:
       muMuTauTauAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
