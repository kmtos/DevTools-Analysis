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

logger = logging.getLogger("ElectronAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class ElectronAnalysis(AnalysisBase):
    '''
    all electrons
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','eTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ETree')
        super(ElectronAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup analysis tree

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('e')
        self.addDetailedElectron('e')

        # met
        self.addMet('met')

    ####################################################
    ### override analyze to store after every lepton ###
    ####################################################
    def perRowAction(self):
        '''Per row action, can be overridden'''
        for elec in self.electrons:
            cands = {
                'e': elec,
                'met': self.pfmet,
                'w': MetCompositeCandidate(self.pfmet,elec),
                'event':self.event,
            }
            if elec.pt()<10: continue
            self.tree.fill(cands,allowDuplicates=True)

        self.eventsStored += 1


    #########################
    ### detailed electron ###
    #########################
    def addDetailedElectron(self,label):
        '''Add detailed electron variables'''
        self.addCandVar(label,'cutBasedVeto','cutBasedVeto','I')
        self.addCandVar(label,'cutBasedLoose','cutBasedLoose','I')
        self.addCandVar(label,'cutBasedMedium','cutBasedMedium','I')
        self.addCandVar(label,'cutBasedTight','cutBasedTight','I')
        self.addCandVar(label,'cutBasedHLTPreselection','cutBasedHLTPreselection','I')
        self.addCandVar(label,'wwLoose','wwLoose','I')
        #self.addCandVar(label,'heepV70','heepV70','I')
        self.addCandVar(label,'mvaValues','mvaValues','F')
        self.addCandVar(label,'mvaCategories','mvaCategories','I')
        self.addCandVar(label,'mvaWP80','mvaWP80','I')
        self.addCandVar(label,'mvaWP90','mvaWP90','I')
        self.addCandVar(label,'susyMVA','susyMVA','F')
        self.addCandVar(label,'isSUSYMVAPreselection','isSUSYMVAPreselection','I')
        self.addCandVar(label,'isSUSYTight','isSUSYTight','I')
        self.addCandVar(label,'isSUSYVLoose','isSUSYVLoose','I')
        self.addCandVar(label,'miniIsolation','miniIsolation','F')
        self.addCandVar(label,'miniIsolationCharged','miniIsolationCharged','F')
        self.addCandVar(label,'miniIsolationNeutral','miniIsolationNeutral','F')
        self.addCandVar(label,'miniIsolationPhoton','miniIsolationPhoton','F')
        self.addCandVar(label,'miniIsolationPileup','miniIsolationPileup','F')
        self.tree.add(lambda cands: cands[label].dr03EcalRecHitSumEt()/cands[label].pt(), '{0}_ecalRelIso'.format(label), 'F')
        self.tree.add(lambda cands: cands[label].dr03HcalTowerSumEt()/cands[label].pt(), '{0}_hcalRelIso'.format(label), 'F')
        self.tree.add(lambda cands: cands[label].dr03TkSumPt()/cands[label].pt(), '{0}_trackRelIso'.format(label), 'F')
        self.addCandVar(label,'superClusterEta','superClusterEta','F')
        self.addCandVar(label,'sigmaIetaIeta','sigmaIetaIeta','F')
        self.addCandVar(label,'hcalOverEcal','hcalOverEcal','F')
        self.addCandVar(label,'deltaEtaSuperClusterTrackAtVtx','deltaEtaSuperClusterTrackAtVtx','F')
        self.addCandVar(label,'deltaPhiSuperClusterTrackAtVtx','deltaPhiSuperClusterTrackAtVtx','F')
        self.addCandVar(label,'passConversionVeto','passConversionVeto','I')
        self.addCandVar(label,'missingHits','missingHits','I')
        self.addCandVar(label,'ecalEnergy','ecalEnergy','F')
        self.addCandVar(label,'eSuperClusterOverP','eSuperClusterOverP','F')
        self.tree.add(lambda cands: abs(1.-cands[label].eSuperClusterOverP())*1./cands[label].ecalEnergy(), '{0}_oneOverEMinusOneOverP'.format(label), 'F')
        self.addCandVar(label,'relPFIsoRhoR03','relPFIsoRhoR03','F')
        self.addCandVar(label,'jetPtRatio','jetPtRatio','F')
        self.addCandVar(label,'jetPtRel','jetPtRel','F')
        self.addCandVar(label,'jetNumberOfChargedDaughters','jetNumberOfChargedDaughters','F')
        self.addCandVar(label,'jetBtagCSV','jetBtagCSV','F')


    def passMVATrigPre(self,cand):
        pt = cand.pt()
        sceta = cand.superClusterEta()
        sigmaIEtaIEta = cand.sigmaIetaIeta()
        hcalOverEcal = cand.hcalOverEcal()
        ecalRelIso = cand.dr03EcalRecHitSumEt()/pt
        hcalRelIso = cand.dr03HcalTowerSumEt()/pt
        trackRelIso = cand.dr03TkSumPt()/pt
        dEtaSC = cand.deltaEtaSuperClusterTrackAtVtx()
        dPhiSC = cand.deltaPhiSuperClusterTrackAtVtx()
        relIsoRho = cand.relPFIsoRhoR03()
        passConversion = cand.passConversionVeto()
        dxy = cand.dB2D()
        dz = cand.dz()
        ecalEnergy = cand.ecalEnergy()
        eSuperClusterOverP = cand.eSuperClusterOverP()
        ooEmooP = abs((1.-eSuperClusterOverP)*1./ecalEnergy)
        passMVATrigPre = True
        if pt<15: passMVATrigPre = False
        if sceta<1.479:
            if sigmaIEtaIEta>0.012: passMVATrigPre = False
            if hcalOverEcal>0.09:   passMVATrigPre = False
            if ecalRelIso>0.37:     passMVATrigPre = False
            if hcalRelIso>0.25:     passMVATrigPre = False
            if trackRelIso>0.18:    passMVATrigPre = False
            if abs(dEtaSC)>0.0095:  passMVATrigPre = False
            if abs(dPhiSC)>0.065:   passMVATrigPre = False
        else:
            if sigmaIEtaIEta>0.033: passMVATrigPre = False
            if hcalOverEcal>0.09:   passMVATrigPre = False
            if ecalRelIso>0.45:     passMVATrigPre = False
            if hcalRelIso>0.28:     passMVATrigPre = False
            if trackRelIso>0.18:    passMVATrigPre = False
        return passMVATrigPre



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='eTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    eAnalysis = ElectronAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ETree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       eAnalysis.analyze()
       eAnalysis.finish()
    except KeyboardInterrupt:
       eAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
