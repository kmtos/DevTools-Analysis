#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *
import KinematicFitter

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MuMuAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class MuMuAnalysis(AnalysisBase):
    '''
    MuMu analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','muMuTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MuMuTree')
        super(MuMuAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passMuon)), 'numMuons', 'I')

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

        # mm leptons
        self.addDiLepton('mm')
        self.addCompositeMet('mmmet')
        self.addLepton('m1')
        self.addDetailedMuon('m1')
        self.addLeptonMet('m1met')
        self.addLepton('m2')
        self.addDetailedMuon('m2')
        self.addLeptonMet('m2met')

        # met
        self.addMet('met')

    ############################
    ### Additional functions ###
    ############################
    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        self.addCandVar(label,'matches_IsoMu24','matches_IsoMu24','I')
        self.addCandVar(label,'matches_IsoTkMu24','matches_IsoTkMu24','I')
        super(MuMuAnalysis,self).addDetailedMuon(label)

    def passMuon(self,cand):
        if cand.pt()<3: return False
        if not cand.isPFMuon(): return False
        if not (cand.isGlobalMuon() or cand.isTrackerMuon()): return False
        if abs(cand.dxy())>=0.2: return False
        if abs(cand.dz())>=0.5: return False
        return True


    #########################
    ### select candidates ###
    #########################
    def selectCandidates(self):
        # reset kinfit
        self._kinfit = None

        candidate = {
            'm1' : None,
            'm1met' : None,
            'm2' : None,
            'm2met' : None,
            'mm' : None,
            'mmmet' : None,
            'met': self.pfmet,
        }

        # get leptons
        muons = [m for m in self.muons if self.passMuon(m)]
        if len(muons)<2:
            return candidate


        # get the candidates
        cand = []
        mmDeltaR = 999
        for pair in itertools.permutations(muons,2):
            # trigger match
            matchTrigger = pair[0].matches_IsoMu24() or pair[0].matches_IsoTkMu24()
            if not matchTrigger: continue
            # charge OS
            if pair[0].charge()==pair[1].charge(): continue
            # require lead m pt>25
            if pair[0].pt()<25: continue
            # make composites
            mm = DiCandidate(pair[0],pair[1])
            # choose best
            if not cand: cand = pair
            better = True
            if mm.deltaR()>mmDeltaR:
                better = False
            if better:
                cand = pair
                mmDeltaR = mm.deltaR()

        if not cand:
            return candidate

        m1 = cand[0] if cand[0].pt()>cand[1].pt() else cand[1]
        m2 = cand[1] if cand[0].pt()>cand[1].pt() else cand[0]

        mm = DiCandidate(m1,m2)
        if mm.M()>30:
            return candidate
        #if mm.deltaR()>1.5: return candidate

        candidate['m1'] = m1
        candidate['m1met'] = MetCompositeCandidate(self.pfmet,m1)
        candidate['m2'] = m2
        candidate['m2met'] = MetCompositeCandidate(self.pfmet,m2)
        candidate['mm'] = DiCandidate(m1,m2)
        candidate['mmmet'] = MetCompositeCandidate(self.pfmet,m1,m2)

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
        result = self.checkTrigger(*datasets,**triggerNames)
        if not result: self.report_failure('fails trigger')
        return result

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['m1','m2']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)







def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('SingleMuon',version='80XMuMu',n=10), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='muMuTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    muMuAnalysis = MuMuAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MuMuTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       muMuAnalysis.analyze()
       muMuAnalysis.finish()
    except KeyboardInterrupt:
       muMuAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
