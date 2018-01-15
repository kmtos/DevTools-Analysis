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
        self.cutTree.add(self.metFilter,'metFilter')
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
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
            #self.tree.add(lambda cands: self.event.Mu45_eta2p1Pass(), 'pass_Mu45_eta2p1', 'I')
            #self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.Ele25_eta2p1_WPTight_GsfPass(), 'pass_Ele25_eta2p1_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')


        # z leptons
        self.addDiLepton('z')
        #self.addDiCandVar('z','z1','z2','mass_uncorrected','mass','F',uncorrected=True)
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z1'],cands['z2']), 'z_zeppenfeld','F')
        self.addLepton('z1')
        self.tree.add(lambda cands: self.passMedium(cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z1'])[0], 'z1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z1'])[0], 'z1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z1'])[0], 'z1_tightScale', 'F')
        self.addLepton('z2')
        self.tree.add(lambda cands: self.passMedium(cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z2'])[0], 'z2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z2'])[0], 'z2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z2'])[0], 'z2_tightScale', 'F')

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
        medLeps = self.getPassingCands('Medium',self.electrons,self.muons)
        if len(medLeps)<2: return candidate # need at least 2 leptons

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for zpair in itertools.combinations(medLeps,2):
            if zpair[0].collName!=zpair[1].collName: continue # SF
            if zpair[0].charge()==zpair[1].charge(): continue # OS
            z = [zpair[0],zpair[1]] if zpair[0].pt()>zpair[1].pt() else [zpair[1],zpair[0]]
            if not bestZ: bestZ = z
            if z[0].pt()>bestZ[0].pt():
                bestZ = z
            if z[0].pt()==bestZ[0].pt() and z[1].pt()>bestZ[1].pt():
                bestZ = z
            #z = DiCandidate(*zpair)
            #massdiff = abs(z.M()-ZMASS)
            #if massdiff<bestMassdiff:
            #    bestZ = zpair
            #    bestMassdiff = massdiff

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])

        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps+[z[0],z[1]],0.4)

        return candidate

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
        return len(self.getPassingCands('Loose',self.electrons,self.muons))>=2

    def trigger(self,cands):
        # apply trigger in data and mc
        #if self.event.isData()<0.5: return True
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
                    #'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL',
                    #'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL',
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'SingleMuon'     : [
                    #'IsoMu22',
                    #'IsoTkMu22',
                    'IsoMu24',
                    'IsoTkMu24',
                    #'Mu45_eta2p1',
                    #'Mu50',
                ],
                'SingleElectron' : [
                    'Ele27_WPTight_Gsf',
                ],
            }


        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        if not cands['z1']: return False
        if cands['z1'].__class__.__name__=='Electron':
            datasets = [
                'DoubleEG', 
                'SingleElectron',
            ]
        else:
            datasets = [
                'DoubleMuon', 
                'SingleMuon',
            ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['z1','z2']]
        if candList[0].__class__.__name__=='Electron':
            triggerList = ['Ele23_WPLoose','Ele17_Ele12'] if self.version=='76X' else ['Ele27Tight','Ele23Ele12']
        else:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Mu17_Mu8'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Mu17Mu8']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList)




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
