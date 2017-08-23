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

logger = logging.getLogger("LLGAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class LLGAnalysis(AnalysisBase):
    '''
    LLG analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','llgTree.root')
        outputTreeName = kwargs.pop('outputTreeName','LLGTree')
        # setup a preselection
        self.preselection = '1'
        super(LLGAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)


        # setup cut tree
        #self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree
        self.tree.add(self.getChannelString, 'channel', ['C',4])

        # trigger
        self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')


        # z leptons
        self.addDiLepton('ll')
        self.addLepton('l1',doId=True,doScales=True)
        self.addLepton('l2',doId=True,doScales=True)

        self.addDiCandidate('lg1')
        self.addDiCandidate('lg2')

        self.addComposite('llg')
        self.addPhoton('g',doId=True,doScales=True)

        # met
        self.addMet('met')

    ############################
    ### select LLG candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'l1' : None,
            'l2' : None,
            'g' : None,
            'll' : None,
            'lg1' : None,
            'lg2' : None,
            'llg' : None,
            'met': self.pfmet,
        }

        # get leptons
        #leps = self.getPassingCands('Loose',self.electrons,self.muons)
        leps = self.electrons+self.muons
        if len(leps)<2: return candidate # need at least 2 leptons
        #phos = self.getPassingCands('PhotonPreselectionNoElectronVeto',self.photons)
        phos = self.photons
        if len(phos)<1: return candidate

        # get invariant masses
        best = ()
        for ll in itertools.combinations(leps,2):
            l1, l2 = ll
            if l1.collName!=l2.collName: continue # SF
            if l1.charge()==l2.charge(): continue # OS
            if l1.pt()<l2.pt(): continue
            if best and l1.pt()<best[0].pt():
                continue
            elif best and l1.pt()==best[0].pt() and l2.pt()<best[2].pt():
                continue
            g = phos[0]
            best = (l1, l2, g)

        if not best: return candidate # need a z candidate


        l1,l2,g = best
        ll = DiCandidate(l1,l2)
        lg1 = DiCandidate(l1,g)
        lg2 = DiCandidate(l2,g)
        llg = CompositeCandidate(l1,l2,g)

        candidate['l1'] = l1
        candidate['l2'] = l2
        candidate['ll'] = ll
        candidate['lg1'] = lg1
        candidate['lg2'] = lg2
        candidate['g'] = g
        candidate['llg'] = llg

        return candidate

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['l1','l2','g']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def twoLoose(self,cands):
        return len(self.getPassingCands('Loose',self.electrons))>=2

    def trigger(self,cands):
        # apply trigger in data and mc
        #if self.event.isData()<0.5: return True
        triggerNames = {
            'DoubleMuon'     : [
                'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
            ],
            'DoubleEG'     : [
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ'
            ],
        }


        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon', 
            'DoubleEG', 
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['l1','l2']]
        triggerList = ['Ele23Ele12'] if isinstance(candList[0],Electron) else ['Mu17Mu8']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList)




def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy',version='80XPhoton'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='llgTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    llgAnalysis = LLGAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='LLGTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       llgAnalysis.analyze()
       llgAnalysis.finish()
    except KeyboardInterrupt:
       llgAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
