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

sys.argv.append('-b')
import ROOT
sys.argv.pop()

logger = logging.getLogger("SingleJetAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class SingleJetAnalysis(AnalysisBase):
    '''
    Select a single lepton to perform a dijet control fake rate
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','singleJetTree.root')
        outputTreeName = kwargs.pop('outputTreeName','SingleJetTree')
        super(SingleJetAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup analysis tree

        # trigger
        self.tree.add(lambda cands: self.event.Photon175Pass(), 'pass_Photon175', 'I')

        # lead jet
        self.addJet('jet')

        # lepton
        self.addPhoton('pho')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            'jet' : None,
            'pho' : Candidate(None),
            'met' : self.pfmet
        }

        # add jet
        jets = self.getCands(self.jets, lambda cand: cand.isLoose()>0.5)
        if len(jets)<1: return candidate
        if jets[0].pt()<175: return candidate

        # get photons
        gs = self.getPassingCands('PhotonPreselection',self.photons)
        cleanGs = self.cleanCands(gs,jets,0.5)

        candidate['jet'] = jets[0]
        if len(cleanGs): candidate['pho'] = cleanGs[0]

        return candidate


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('qcd',version='80XQCD'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='singleJetTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    singleJetAnalysis = SingleJetAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='SingleJetTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       singleJetAnalysis.analyze()
       singleJetAnalysis.finish()
    except KeyboardInterrupt:
       singleJetAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
