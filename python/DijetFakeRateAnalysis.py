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

logger = logging.getLogger("DijetFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class DijetFakeRateAnalysis(AnalysisBase):
    '''
    Select a single lepton to perform a dijet control fake rate
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','dijetFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','DijetFakeRateTree')
        super(DijetFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.vetoSecond,'vetoSecond')
        self.cutTree.add(self.metVeto,'metVeto')
        self.cutTree.add(self.mtVeto,'mtVeto')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',2])

        # trigger
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVLPass(), 'pass_Mu8_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Ele23_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele23_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')
        self.tree.add(self.triggerPrescale, 'triggerPrescale', 'F')

        # lead jet
        self.addJet('leadJet')
        self.addDiCandidate('lj')

        # lepton
        self.addLeptonMet('w')
        self.addLepton('l1')
        self.tree.add(lambda cands: self.passMedium(cands['l1']), 'l1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['l1']), 'l1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['l1'])[0], 'l1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['l1'])[0], 'l1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['l1'])[0], 'l1_tightScale', 'F')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            'l1' : None,
            'w' : None,
            'leadJet' : None,
            'lj': None,
            'met': self.pfmet
        }

        # get leptons
        leps = self.getPassingCands('Loose',self.electrons,self.muons)
        if len(leps)<1: return candidate # need at least 1 lepton

        # choose highest pt
        bestCand = 0
        bestPt = 0
        for l in leps:
           if l.pt()>bestPt:
               bestPt = l.pt()
               bestCand = l

        # add jet
        jets = self.getCands(self.jets, lambda cand: cand.isLoose()>0.5)
        if len(jets)<1: return candidate # need a recoil jet

        # choose highest pt
        bestJet = 0
        bestPt = 0
        for j in jets:
           if j.pt()>bestPt:
               bestPt = j.pt()
               bestJet = j

        candidate['l1'] = bestCand
        candidate['w'] = MetCompositeCandidate(self.pfmet,bestCand)
        candidate['leadJet'] = bestJet
        candidate['lj'] = DiCandidate(bestCand,bestJet)

        return candidate


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['l1']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def vetoSecond(self,cands):
        return len(self.getPassingCands('Loose',self.electrons,self.muons))==1

    def metVeto(self,cands):
        return cands['met'].et()<25

    def mtVeto(self,cands):
        if not cands['w']: return False
        return cands['w'].Mt()<25

    def trigger(self,cands):
        # accept MC, check trigger for data
        isData = self.event.isData()>0.5
        if not cands['l1']: return False
        triggerNames = {
            'DoubleMuon'     : [
                ['Mu8_TrkIsoVVL', 0],
                ['Mu17_TrkIsoVVL', 25],
            ],
            'DoubleEG'       : [
                ['Ele12_CaloIdL_TrackIdL_IsoVL', 0],
                ['Ele23_CaloIdL_TrackIdL_IsoVL', 30],
            ],
        }
        # here we need to accept only a specific trigger after a certain pt threshold
        pt = cands['l1'].pt()
        dataset = 'DoubleEG' if isinstance(cands['l1'],Electron) else 'DoubleMuon'
        # accept the event only if it is triggered in the current dataset
        reject = True if isData else False
        if dataset in self.fileNames[0]: reject = False
        # now pick the appropriate trigger for the pt
        theTrigger = ''
        for trig, ptThresh in triggerNames[dataset]:
            if pt < ptThresh: break
            theTrigger = trig
        # and check if we pass
        passTrigger = getattr(self.event,'{0}Pass'.format(theTrigger))()
        if passTrigger>0.5: return False if reject else True
        return False

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands['l1']]
        # select via pt and flavor
        pt = cands['l1'].pt()
        if isinstance(cands['l1'],Electron):
            #triggerList = ['Ele17_Ele12Leg2'] if self.version=='76X' else ['Ele23Ele12Leg2']
            if pt<25:
                triggerList = ['Ele17_Ele12Leg2'] if self.version=='76X' else ['Ele23Ele12Leg2']
            else:
                triggerList = ['Ele17_Ele12Leg1'] if self.version=='76X' else ['Ele23Ele12Leg1']
        else:
            #triggerList = ['Mu17_Mu8Leg2'] if self.version=='76X' else ['Mu17Mu8Leg2']
            if pt<25:
                triggerList = ['Mu17_Mu8Leg2'] if self.version=='76X' else ['Mu17Mu8Leg2']
            else:
                triggerList = ['Mu17_Mu8Leg1'] if self.version=='76X' else ['Mu17Mu8Leg1']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList)

    def triggerPrescale(self,cands):
        # select via pt and flavor
        pt = cands['l1'].pt()
        if isinstance(cands['l1'],Electron):
            #trigger = 'Ele17_Ele12Leg2'
            if pt<25:
                trigger = 'Ele23_Ele12Leg2'
            else:
                trigger = 'Ele23_Ele12Leg1'
        else:
            #trigger = 'Mu17_Mu8Leg2'
            if pt<25:
                trigger = 'Mu17_Mu8Leg2'
            else:
                trigger = 'Mu17_Mu8Leg1'
        return self.triggerPrescales.getPrescale(trigger)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('data'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='dijetFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dijetFakeRateAnalysis = DijetFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='DijetFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       dijetFakeRateAnalysis.analyze()
       dijetFakeRateAnalysis.finish()
    except KeyboardInterrupt:
       dijetFakeRateAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
