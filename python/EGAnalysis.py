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

logger = logging.getLogger("EGAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class EGAnalysis(AnalysisBase):
    '''
    WZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','egTree.root')
        outputTreeName = kwargs.pop('outputTreeName','EGTree')
        super(EGAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')

        # trigger
        #self.addPhotonTriggers()
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # electron photon
        self.addDiCandidate('eg')

        # z leptons
        self.addLepton('e',doId=True,doScales=True)
        self.addPhoton('g',doId=True)

        # met
        self.addMet('met')

    ############################
    ### select WZ candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'e' : None,
            'g' : None,
            'eg': None,
            'met': self.pfmet,
            'cleanJets': [],
            #'leadJet': Candidate(None),
            #'subleadJet': Candidate(None),
            #'dijet': DiCandidate(Candidate(None),Candidate(None)),
        }

        gs = self.getPassingCands('PhotonPreselectionNoElectronVeto',self.photons)
        es = self.getPassingCands('Medium',self.electrons)

        if not es or not gs: return candidate

        best = ()
        for e in es:
            if e.pt()<30: continue
            for g in gs:
                if deltaR(e.eta(),e.phi(),g.eta(),g.phi())<0.4: continue
                if not best: best = (e,g)
                if e.pt()>best[0].pt():
                    best = (e,g)
                elif g.pt()>best[1].pt():
                    best = (e,g)

        if not best: return candidate
        e,g = best

        candidate['e'] = e
        candidate['g'] = g
        candidate['eg'] = DiCandidate(e,g)

        goodGs = self.getPassingCands('Photon',self.photons)
        candidate['cleanJets'] = self.cleanCands(self.jets,es+goodGs+[g],0.4)

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
            'SingleElectron' : [
                'Ele27_WPTight_Gsf',
            ],
        }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'SingleElectron',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands['e']]
        triggerList = ['Ele27Tight']
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
    parser.add_argument('--outputFile', type=str, default='egTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    egAnalysis = EGAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='EGTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       egAnalysis.analyze()
       egAnalysis.finish()
    except KeyboardInterrupt:
       egAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
