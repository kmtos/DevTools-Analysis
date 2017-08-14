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

logger = logging.getLogger("MMGAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class MMGAnalysis(AnalysisBase):
    '''
    MMG analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','mmgTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MMGTree')
        # setup a preselection
        self.preselection = 'muons_count>1'
        super(MMGAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)


        # setup cut tree
        self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # trigger
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')


        # z leptons
        self.addDiLepton('z')
        self.addLepton('z1',doId=True,doScales=True)
        self.addLepton('z2',doId=True,doScales=True)

        self.addComposite('zg')
        self.addPhoton('g',doId=True,doScales=True)

        # met
        self.addMet('met')

    ############################
    ### select MMG candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            'g' : None,
            'z' : None,
            'zg' : None,
            'met': self.pfmet,
        }

        # get leptons
        medLeps = self.getPassingCands('Medium',self.muons)
        if len(medLeps)<2: return candidate # need at least 2 leptons
        #phos = self.getPassingCands('PhotonPreselectionNoElectronVeto',self.photons)
        phos = self.photons

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for mm in itertools.combinations(medLeps,2):
            m1, m2 = mm
            if m1.collName!=m2.collName: continue # SF
            if m1.charge()==m2.charge(): continue # OS
            for g in phos:
                if deltaR(m1.eta(),m1.phi(),g.eta(),g.phi())<0.8:
                    if m2.pt()<20: continue
                elif deltaR(m2.eta(),m2.phi(),g.eta(),g.phi())<0.8:
                    if m1.pt()<20: continue
                else:
                    continue
                zg = CompositeCandidate(m1,m2,g)
                massdiff = abs(zg.M()-ZMASS)
                if massdiff<bestMassdiff:
                    bestZ = (m1,m2,g) if m1.pt()>m2.pt() else (m2,m1,g)
                    bestMassdiff = massdiff

        if not bestZ: return candidate # need a z candidate

        m1,m2,g = bestZ
        candidate['z1'] = m1
        candidate['z2'] = m2
        candidate['z'] = DiCandidate(m1,m2)
        candidate['g'] = g
        candidate['zg'] = CompositeCandidate(m1,m2,g)

        return candidate

    ###########################
    ### analysis selections ###
    ###########################
    def twoLoose(self,cands):
        return len(self.getPassingCands('Loose',self.muons))>=2

    def trigger(self,cands):
        # apply trigger in data and mc
        #if self.event.isData()<0.5: return True
        triggerNames = {
            'DoubleMuon'     : [
                'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
            ],
        }


        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon', 
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['z1','z2']]
        triggerList = ['Mu17_Mu8'] if self.version=='76X' else ['Mu17Mu8']
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
    parser.add_argument('--outputFile', type=str, default='mmgTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    mmgAnalysis = MMGAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MMGTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       mmgAnalysis.analyze()
       mmgAnalysis.finish()
    except KeyboardInterrupt:
       mmgAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
