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

logger = logging.getLogger("TwoPhotonAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class TwoPhotonAnalysis(AnalysisBase):
    '''
    WZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','twoPhotonTree.root')
        outputTreeName = kwargs.pop('outputTreeName','TwoPhotonTree')
        super(TwoPhotonAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')

        # trigger
        self.addPhotonTriggers()
        #self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        #self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        #self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # z leptons
        self.addDiCandidate('gg')
        self.addPhoton('g1',doId=True)
        self.addPhoton('g2',doId=True)

        # met
        self.addMet('met')

    ############################
    ### select WZ candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'g1' : None,
            'g2' : None,
            'gg' : None,
            'met': self.pfmet,
            'cleanJets': [],
            #'leadJet': Candidate(None),
            #'subleadJet': Candidate(None),
            #'dijet': DiCandidate(Candidate(None),Candidate(None)),
        }

        gs = self.getPassingCands('PhotonPreselection',self.photons)

        best = ()

        for gCand in itertools.permutations(gs,2):
            g1, g2 = gCand
            if g1.pt()<g2.pt(): continue
            if g1.pt()<32: continue
            if g2.pt()<20: continue
            if not best: best = gCand
            if g1.pt()>best[0].pt():
                best = gCand
            elif g1.pt()==best[0].pt() and g2.pt()>best[1].pt():
                best = gCand
            

        if not best: return candidate

        g1, g2 = best
        candidate['g1'] = g1
        candidate['g2'] = g2
        candidate['gg'] = DiCandidate(g1,g2)

        goodGs = self.getPassingCands('Photon',self.photons)
        candidate['cleanJets'] = self.cleanCands(self.jets,goodGs+[g1,g2],0.4)

        return candidate


    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        # accept MC, check trigger for data
        #if self.event.isData()<0.5: return True
        # use trigger for MC and data
        isData = self.event.isData()>0.5
        triggerNames = {
            'DoubleEG'       : [
                'DoublePhoton60',
                'Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90',
            ],
            'SinglePhoton' : [
                'Photon175',
            ],
        }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleEG', 
            'SinglePhoton',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    #def triggerEfficiencyMC(self,cands):
    #    return self.triggerEfficiency(cands,mode='mc')

    #def triggerEfficiencyData(self,cands):
    #    return self.triggerEfficiency(cands,mode='data')

    #def triggerEfficiency(self,cands,mode='ratio'):
    #    candList = [cands[c] for c in ['g1','g2','g3']]
    #    triggerList = []
    #    if mode=='data':
    #        return self.triggerScales.getDataEfficiency(triggerList,candList)
    #    elif mode=='mc':
    #        return self.triggerScales.getMCEfficiency(triggerList,candList)
    #    elif mode=='ratio':
    #        return self.triggerScales.getRatio(triggerList,candList)








def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy',version='80XPhoton'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='twoPhotonTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    twoPhotonAnalysis = TwoPhotonAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='TwoPhotonTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       twoPhotonAnalysis.analyze()
       twoPhotonAnalysis.finish()
    except KeyboardInterrupt:
       twoPhotonAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
