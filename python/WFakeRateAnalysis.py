#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import sys
import itertools
import operator

sys.argv.append('-b')
import ROOT
sys.argv.pop()

logger = logging.getLogger("WFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class WFakeRateAnalysis(AnalysisBase):
    '''
    Select a muon + lep for fake rate in a W control region
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','WFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','WFakeRateTree')
        super(WFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu22Pass(), 'pass_IsoMu22', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu22Pass(), 'pass_IsoTkMu22', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # lepton
        # mu tag
        self.addLeptonMet('wm')
        self.addLepton('m')
        self.tree.add(lambda cands: self.tightScale(cands['m']), 'm_tightScale', 'F')

        # lep probe
        self.addLeptonMet('wl')
        self.addLepton('l')
        self.tree.add(lambda cands: self.passMedium(cands['l']), 'l_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['l']), 'l_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['l']), 'l_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['l']), 'l_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['l']), 'l_tightScale', 'F')

        # dilepton combination
        self.addDiLepton('z')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            'm' : None,
            'l' : None,
            'z' : None,
            'wm': None,
            'wl': None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        muons = self.getCands(self.muons,self.passTight)
        if len(muons)!=1: return candidate # need 1 tight muon
        if len(leps)!=2: return candidate # need 1 additional lep

        # get invariant masses
        bestL = ()
        bestPt = 0
        for z2 in leps:
            zpair = (muons[0], z2)
            if zpair[0].pt()<25: continue
            if zpair[1].pt()<20: continue
            z = DiCandidate(*zpair)
            if z.deltaR()<0.5: continue
            pt = zpair[1].pt()
            if pt>bestPt:
                bestL = zpair
                bestPt = pt

        if not bestL: return candidate # need a z candidate



        z = bestL
        candidate['m'] = z[0]
        candidate['l'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['wm'] = MetCompositeCandidate(self.pfmet,z[0])
        candidate['wl'] = MetCompositeCandidate(self.pfmet,z[1])


        return candidate


    ##################
    ### lepton IDs ###
    ##################
    def passLoose(self,cand):
        return passHppLoose(cand)

    def passMedium(self,cand):
        return passHppMedium(cand)

    def passTight(self,cand):
        return passHppTight(cand)

    def looseScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDLooseIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedVeto',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDVeto',cand)
        else:                           return 1.

    def mediumScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedMedium',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDMedium',cand)
        else:                           return 1.

    def tightScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedTight',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDTight',cand)
        else:                           return 1.

    def getPassingCands(self,mode):
        if mode=='Loose':
            passMode = self.passLoose
        elif mode=='Medium':
            passMode = self.passMedium
        elif mode=='Tight':
            passMode = self.passTight
        else:
            return []
        cands = []
        for coll in [self.muons,self.electrons]:
            cands += self.getCands(coll,passMode)
        return cands


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = 'm' + self.getCollectionString(cands['l'])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
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
                    'IsoMu22',
                    'IsoTkMu22',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'SingleMuon',
        ]
        # reject triggers if they are in another dataset
        # looks for the dataset name in the filename
        # for MC it accepts all
        reject = True if self.event.isData()>0.5 else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0]: reject = False
            for trigger in triggerNames[dataset]:
                var = '{0}Pass'.format(trigger)
                passTrigger = getattr(self.event,var)()
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0]: break
        return False

    def triggerEfficiency(self,cands):
        candList = [cands['m']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu22ORIsoTkMu22']
        return self.triggerScales.getDataEfficiency(triggerList,candList)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('SingleMuon'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='wFakeRateTree.root', help='Output file')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = WFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='WFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
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