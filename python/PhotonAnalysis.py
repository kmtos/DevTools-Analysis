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

logger = logging.getLogger("PhotonAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class PhotonAnalysis(AnalysisBase):
    '''
    all electrons
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','gTree.root')
        outputTreeName = kwargs.pop('outputTreeName','GTree')
        super(PhotonAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup analysis tree

        self.addPhoton('g',doId=True)
        self.addDetailedPhoton('g')


    ####################################################
    ### override analyze to store after every lepton ###
    ####################################################
    def perRowAction(self):
        '''Per row action, can be overridden'''
        for pho in self.photons:
            cands = {
                'g': pho,
                'event':self.event,
            }
            if pho.pt()<10: continue
            self.tree.fill(cands,allowDuplicates=True)

        self.eventsStored += 1


    #######################
    ### detailed photon ###
    #######################
    def addDetailedPhoton(self,label):
        '''Add detailed photon variables'''
        self.addCandVar(label,'superClusterEta','superClusterEta','F')
        self.addCandVar(label,'superClusterPhi','superClusterPhi','F')
        self.addCandVar(label,'superClusterEnergy','superClusterEnergy','F')
        self.addCandVar(label,'superClusterRawEnergy','superClusterRawEnergy','F')
        self.addCandVar(label,'superClusterEtaWidth','superClusterEtaWidth','F')
        self.addCandVar(label,'superClusterPhiWidth','superClusterPhiWidth','F')
        self.addCandVar(label,'hadronicOverEM','hadronicOverEM','F')
        self.addCandVar(label,'hadronicDepth1OverEm','hadronicDepth1OverEm','F')
        self.addCandVar(label,'hadronicDepth2OverEm','hadronicDepth2OverEm','F')
        self.addCandVar(label,'sigmaIEtaIEta','sigmaIEtaIEta','F')
        self.addCandVar(label,'e1x5','e1x5','F')
        self.addCandVar(label,'e2x5','e2x5','F')
        self.addCandVar(label,'e3x3','e3x3','F')
        self.addCandVar(label,'e5x5','e5x5','F')
        self.addCandVar(label,'maxEnergyXtal','maxEnergyXtal','F')
        self.addCandVar(label,'r1x5','r2x5','F')
        self.addCandVar(label,'r2x5','r2x5','F')
        self.addCandVar(label,'gammaDR030','gammaDR030','F')
        self.addCandVar(label,'phoWorstChargedIsolationWithConeVeto','phoWorstChargedIsolationWithConeVeto','F')
        self.addCandVar(label,'phoESEffSigmaRR','phoESEffSigmaRR','F')
        self.addCandVar(label,'phoFull5x5E1x3','phoFull5x5E1x3','F')
        self.addCandVar(label,'phoFull5x5E2x2','phoFull5x5E2x2','F')
        self.addCandVar(label,'phoFull5x5E2x5Max','phoFull5x5E2x5Max','F')
        self.addCandVar(label,'phoFull5x5E5x5','phoFull5x5E5x5','F')
        self.addCandVar(label,'phoFull5x5SigmaIEtaIEta','phoFull5x5SigmaIEtaIEta','F')
        self.addCandVar(label,'phoFull5x5SigmaIEtaIPhi','phoFull5x5SigmaIEtaIPhi','F')
        self.addCandVar(label,'phoChargedIsolation','phoChargedIsolation','F')
        self.addCandVar(label,'phoNeutralHadronIsolation','phoNeutralHadronIsolation','F')
        self.addCandVar(label,'phoPhotonIsolation','phoPhotonIsolation','F')
        self.addCandVar(label,'effectiveAreaChargedHadrons','effectiveAreaChargedHadrons','F')
        self.addCandVar(label,'effectiveAreaNeutralHadrons','effectiveAreaNeutralHadrons','F')
        self.addCandVar(label,'effectiveAreaPhotons','effectiveAreaPhotons','F')
        self.addCandVar(label,'trackIso','trackIso','F')
        self.tree.add(lambda cands: self.matchToGenPhotons(cands[label]), 'numMatchedPhotons', 'I')

    def matchToGenPhotons(self,pho):
        numMatch = 0
        for gen in self.gen:
            if gen.pdgId() != 22: continue
            if gen.status() != 1: continue
            #if abs(gen.mother_1()) not in [1,2,3,4,5,6,21,22,35,36]: continue # from quarks for gjet mc, 35/36 for signal
            if deltaR(gen.eta(),gen.phi(),pho.eta(),pho.phi())>0.1: continue
            numMatch += 1
        return numMatch





def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('gjet',version='80XPhoton'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='gTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    gAnalysis = PhotonAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='GTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       gAnalysis.analyze()
       gAnalysis.finish()
    except KeyboardInterrupt:
       gAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
