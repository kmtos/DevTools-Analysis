#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import itertools
import operator

import ROOT

logger = logging.getLogger("TauChargeAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class TauChargeAnalysis(AnalysisBase):
    '''
    TauCharge analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','tauChargeTree.root')
        outputTreeName = kwargs.pop('outputTreeName','TauChargeTree')
        super(TauChargeAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('m')
        self.tree.add(lambda cands: self.passMedium(cands['m']), 'm_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['m']), 'm_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['m'])[0], 'm_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['m'])[0], 'm_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['m'])[0], 'm_tightScale', 'F')
        self.addLepton('t')
        self.tree.add(lambda cands: self.passMedium(cands['t']), 't_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['t']), 't_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['t'])[0], 't_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['t'])[0], 't_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['t'])[0], 't_tightScale', 'F')

        # w lepton
        self.addLeptonMet('wm')
        self.addLeptonMet('wt')

        # met
        self.addMet('met')

    ############################
    ### select DY candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'm' : None,
            't' : None,
            'z' : None,
            'wm': None,
            'wt': None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose',self.muons,self.taus)
        if len(leps)<2: return candidate # need at least 2 leptons

        # get invariant masses
        bestZ = ()
        bestst = 0
        for zpair in itertools.permutations(leps,2):
            if zpair[0].collName!='muons': continue
            if zpair[1].collName!='taus': continue
            if zpair[0].pt()<20: continue
            if zpair[1].pt()<20: continue
            z = DiCandidate(*zpair)
            if z.deltaR()<0.02: continue
            st = z.st()
            if st>bestst:
                bestZ = zpair
                bestst = st

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = bestZ
        candidate['m'] = z[0]
        candidate['t'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['wm'] = MetCompositeCandidate(self.pfmet,z[0])
        candidate['wt'] = MetCompositeCandidate(self.pfmet,z[1])


        return candidate


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['m','t']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def twoLoose(self,cands):
        return len(self.getPassingCands('Loose',self.muons,self.taus))>=2

    def trigger(self,cands):
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

    def triggerEfficiency(self,cands):
        candList = [cands['m']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24']
        return self.triggerScales.getDataEfficiency(triggerList,candList)









def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='tauChargeTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = TauChargeAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='TauChargeTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       dyAnalysis.analyze()
       dyAnalysis.finish()
    except KeyboardInterrupt:
       dyAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)

