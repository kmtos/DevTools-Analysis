#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from Candidates import *

import itertools
import operator

import ROOT

logger = logging.getLogger("DYAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class DYAnalysis(AnalysisBase):
    '''
    DY analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','dyTree.root')
        outputTreeName = kwargs.pop('outputTreeName','DYTree')
        # setup a preselection
        self.preselection = '(electrons_count>1 || muons_count>1)'
        super(DYAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)


        # setup cut tree
        self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # alt pileupweights
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,60000), 'pileupWeight_60000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,61000), 'pileupWeight_61000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,62000), 'pileupWeight_62000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,63000), 'pileupWeight_63000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,64000), 'pileupWeight_64000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,65000), 'pileupWeight_65000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,66000), 'pileupWeight_66000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,67000), 'pileupWeight_67000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,68000), 'pileupWeight_68000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,69000), 'pileupWeight_69000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,70000), 'pileupWeight_70000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,71000), 'pileupWeight_71000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,72000), 'pileupWeight_72000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,73000), 'pileupWeight_73000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,74000), 'pileupWeight_74000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,75000), 'pileupWeight_75000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,76000), 'pileupWeight_76000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,77000), 'pileupWeight_77000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,78000), 'pileupWeight_78000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,79000), 'pileupWeight_79000', 'F')
        self.tree.add(lambda cands: self.pileupWeights.alt_weight(self.event,80000), 'pileupWeight_80000', 'F')

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passMedium)), 'numMediumElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passMedium)), 'numMediumMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
            self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu22Pass(), 'pass_IsoMu22', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu22Pass(), 'pass_IsoTkMu22', 'I')
            self.tree.add(lambda cands: self.event.Mu45_eta2p1Pass(), 'pass_Mu45_eta2p1', 'I')
            self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.Ele25_eta2p1_WPTight_GsfPass(), 'pass_Ele25_eta2p1_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_eta2p1_WPLoose_GsfPass(), 'pass_Ele27_eta2p1_WPLoose_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele45_WPLoose_GsfPass(), 'pass_Ele45_WPLoose_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')


        # z leptons
        self.addDiLepton('z')
        #self.addDiCandVar('z','z1','z2','mass_uncorrected','mass','F',uncorrected=True)
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z1'],cands['z2']), 'z_zeppenfeld','F')
        self.addLepton('z1')
        self.tree.add(lambda cands: self.passMedium(cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z1']), 'z1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z1']), 'z1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z1']), 'z1_tightScale', 'F')
        self.addLepton('z2')
        self.tree.add(lambda cands: self.passMedium(cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z2']), 'z2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z2']), 'z2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z2']), 'z2_tightScale', 'F')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select DY candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            'z' : None,
            'met': self.pfmet,
            'cleanJets': [],
        }

        # get leptons
        #leps = self.getPassingCands('Loose')
        medLeps = self.getPassingCands('Medium')
        if len(medLeps)<2: return candidate # need at least 2 leptons

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for zpair in itertools.combinations(medLeps,2):
            if zpair[0].collName!=zpair[1].collName: continue # SF
            if zpair[0].charge()==zpair[1].charge(): continue # OS
            z = DiCandidate(*zpair)
            massdiff = abs(z.M()-ZMASS)
            if massdiff<bestMassdiff:
                bestZ = zpair
                bestMassdiff = massdiff

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])

        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps,0.4)

        return candidate

    #################
    ### lepton id ###
    #################
    def passLoose(self,cand):
        return passHppLoose(cand)

    def passMedium(self,cand):
        return passHppMedium(cand)

    def passTight(self,cand):
        return passHppTight(cand)

    def looseScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDLooseIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedVeto',cand)
        else:                           return 1.

    def mediumScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedMedium',cand)
        else:                           return 1.

    def tightScale(self,cand):
        if isinstance(cand,Muon):       return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif isinstance(cand,Electron): return self.leptonScales.getScale('CutbasedTight',cand)
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
        for coll in [self.electrons,self.muons]:
            cands += self.getCands(coll,passMode)
        return cands

    def numJets(self,cleanJets,mode,pt):
        jetColl = self.getCands(
            cleanJets,
            lambda cand: getattr(cand,mode)()>0.5 and cand.pt()>pt
        )
        return len(jetColl)

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def twoLoose(self,cands):
        return len(self.getPassingCands('Loose'))>=2

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
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
                'SingleElectron' : [
                    'Ele23_WPLoose_Gsf',
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
            }


        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        if not cands['z1']: return False
        if isinstance(cands['z1'],Electron):
            datasets = [
                'DoubleEG', 
                'SingleElectron',
            ]
        else:
            datasets = [
                'DoubleMuon', 
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
        candList = [cands[c] for c in ['z1','z2']]
        if isinstance(candList[0],Electron):
            triggerList = ['Ele23_WPLoose','Ele17_Ele12'] if self.version=='76X' else ['SingleEleSoup','Ele23Ele12']
        else:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Mu17_Mu8'] if self.version=='76X' else ['SingleMuSoup','Mu17Mu8']
        return self.triggerScales.getDataEfficiency(triggerList,candList)





def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('dy'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='dyTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    dyAnalysis = DYAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='DYTree',
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
