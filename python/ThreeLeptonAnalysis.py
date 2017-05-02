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

import ROOT

logger = logging.getLogger("ThreeLeptonAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class ThreeLeptonAnalysis(AnalysisBase):
    '''
    ThreeLepton analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','threeLeptonTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ThreeLeptonTree')
        self.preselection = 'electrons_count+muons_count+taus_count>2'
        super(ThreeLeptonAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        self.addTriggers()
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')

        # leptons
        self.addLepton('e1')
        self.tree.add(lambda cands: self.passMedium(    cands['e1']),        'e1_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['e1']),        'e1_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['e1'])[0],     'e1_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e1'])[1],     'e1_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e1'])[2],     'e1_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e1'])[0],     'e1_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e1'])[1],     'e1_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e1'])[2],     'e1_tightScaleDown',     'F')
        self.addLepton('e2')
        self.tree.add(lambda cands: self.passMedium(    cands['e2']),        'e2_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['e2']),        'e2_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['e2'])[0],     'e2_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e2'])[1],     'e2_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e2'])[2],     'e2_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e2'])[0],     'e2_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e2'])[1],     'e2_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e2'])[2],     'e2_tightScaleDown',     'F')
        self.addLepton('e3')
        self.tree.add(lambda cands: self.passMedium(    cands['e3']),        'e3_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['e3']),        'e3_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['e3'])[0],     'e3_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e3'])[1],     'e3_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['e3'])[2],     'e3_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e3'])[0],     'e3_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e3'])[1],     'e3_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['e3'])[2],     'e3_tightScaleDown',     'F')

        self.addLepton('m1')
        self.tree.add(lambda cands: self.passMedium(    cands['m1']),        'm1_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['m1']),        'm1_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['m1'])[0],     'm1_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m1'])[1],     'm1_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m1'])[2],     'm1_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m1'])[0],     'm1_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m1'])[1],     'm1_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m1'])[2],     'm1_tightScaleDown',     'F')
        self.addLepton('m2')
        self.tree.add(lambda cands: self.passMedium(    cands['m2']),        'm2_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['m2']),        'm2_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['m2'])[0],     'm2_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m2'])[1],     'm2_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m2'])[2],     'm2_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m2'])[0],     'm2_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m2'])[1],     'm2_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m2'])[2],     'm2_tightScaleDown',     'F')
        self.addLepton('m3')
        self.tree.add(lambda cands: self.passMedium(    cands['m3']),        'm3_passMedium',         'I')
        self.tree.add(lambda cands: self.passTight(     cands['m3']),        'm3_passTight',          'I')
        self.tree.add(lambda cands: self.mediumScale(   cands['m3'])[0],     'm3_mediumScale',        'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m3'])[1],     'm3_mediumScaleUp',      'F')
        self.tree.add(lambda cands: self.mediumScale(   cands['m3'])[2],     'm3_mediumScaleDown',    'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m3'])[0],     'm3_tightScale',         'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m3'])[1],     'm3_tightScaleUp',       'F')
        self.tree.add(lambda cands: self.tightScale(    cands['m3'])[2],     'm3_tightScaleDown',     'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['m3'])[0],     'm3_mediumFakeRate',     'F')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select 3l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'e1' : None,
            'e2' : None,
            'e3' : None,
            'm1' : None,
            'm2' : None,
            'm3' : None,
            'met': self.pfmet,
            'cleanJets' : [],
        }

        # get leptons
        leps = self.getPassingCands('Medium',self.electrons,self.muons)
        if len(leps)<3: return candidate # need at least 3 leptons
        elecs = [cand for cand in leps if cand.collName=='electrons']
        muons = [cand for cand in leps if cand.collName=='muons']

        goodElecs = []
        for i,candi in enumerate(elecs):
            keep = True
            for j,candj in enumerate(elecs):
                if j<=i: continue
                dicand = DiCandidate(candi,candj)
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            goodElecs += [candi]

        goodMuons = []
        for i,candi in enumerate(muons):
            keep = True
            for j,candj in enumerate(muons):
                if j<=i: continue
                dicand = DiCandidate(candi,candj)
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            goodMuons += [candi]

        if len(goodElecs+goodMuons)<3: return candidate

        ems = goodElecs+goodMuons
        pts_ems = sorted([cand.pt() for cand in ems])
        if pts_ems[-1]<30.: return candidate
        if len(ems)>1:
            if pts_ems[-2]<20.: return candidate

        e1 = goodElecs[0] if len(goodElecs)>0 else Candidate(None)
        e2 = goodElecs[1] if len(goodElecs)>1 else Candidate(None)
        e3 = goodElecs[2] if len(goodElecs)>2 else Candidate(None)
        m1 = goodMuons[0] if len(goodMuons)>0 else Candidate(None)
        m2 = goodMuons[1] if len(goodMuons)>1 else Candidate(None)
        m3 = goodMuons[2] if len(goodMuons)>2 else Candidate(None)

        candidate['e1'] = e1
        candidate['e2'] = e2
        candidate['e3'] = e3
        candidate['m1'] = m1
        candidate['m2'] = m2
        candidate['m3'] = m3

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,leps,0.4)

        return candidate

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,cands):
        return len(self.getPassingCands('Loose',self.electrons,self.muons))>=3

    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
        if self.version=='76X':
            triggerNames = {
                'DoubleMuon'     : [
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL',
                    'Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
                'SingleElectron' : [
                    'Ele23_WPLoose_Gsf',
                ],
                'Tau' : [
                    'DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg',
                ],
            }
        else:
            triggerNames = {
                'DoubleMuon'     : [
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    'IsoMu22',
                    'IsoTkMu22',
                    'Mu45_eta2p1',
                    'Mu50',
                ],
                'SingleElectron' : [
                    'Ele25_eta2p1_WPTight_Gsf',
                    'Ele27_WPTight_Gsf',
                    'Ele27_eta2p1_WPLoose_Gsf',
                    'Ele45_WPLoose_Gsf',
                ],
                'Tau' : [
                    'DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg',
                ],
            }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'DoubleMuon', 
            'DoubleEG', 
            'MuonEG',
            'SingleMuon',
            'SingleElectron',
            #'Tau',
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
        candList = [cands[c] for c in ['e1','e2','e3','m1','m2','m3']]
        numTaus = [c.collName for c in candList].count('taus')
        triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['SingleMuSoup','SingleEleSoup','Mu17Mu8','Ele23Ele12']
        return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('hpp3l'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='threeLeptonTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    threeLeptonAnalysis = ThreeLeptonAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ThreeLeptonTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       threeLeptonAnalysis.analyze()
       threeLeptonAnalysis.finish()
    except KeyboardInterrupt:
       threeLeptonAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
