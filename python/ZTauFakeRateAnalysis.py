#!/usr/bin/env python
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from Candidates import *
from DevTools.Analyzer.utilities import getTestFiles

import itertools
import operator
import argparse
import logging
import sys

import ROOT

logger = logging.getLogger("ZTauFakeRateAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class ZTauFakeRateAnalysis(AnalysisBase):
    '''
    ZTauFakeRate analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','dyTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ZTauFakeRateTree')
        super(ZTauFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)


        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',4])

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
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLooseNew)), 'numLooseNewTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passMedium)), 'numMediumTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

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

        # fake lepton
        self.addLeptonMet('w')
        self.addLepton('t')
        self.addDetailedTau('t')
        self.addJet('tjet')
        self.tree.add(lambda cands: self.passLoose(cands['t']), 't_passLoose', 'I')
        self.tree.add(lambda cands: self.passMedium(cands['t']), 't_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['t']), 't_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['t'])[0], 't_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['t'])[0], 't_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['t'])[0], 't_tightScale', 'F')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select DY candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1': None,
            'z2': None,
            'z' : None,
            't': None,
            'tjet': Candidate(None),
            'w' : None,
            'met': self.pfmet,
            'cleanJets': [],
        }

        # get leptons
        leps = self.getPassingCands('Medium',self.electrons,self.muons)
        tleps = self.getPassingCands('LooseNew',self.taus)
        if len(leps)!=2: return candidate # need 3 leptons
        if len(tleps)<1: return candidate # need 3 leptons

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for zpair in itertools.permutations(leps,2):
            # z pass medium
            if not self.passMedium(zpair[0]): continue
            if not self.passMedium(zpair[1]): continue
            if zpair[0].collName!=zpair[1].collName: continue # SF
            if zpair[0].charge()==zpair[1].charge(): continue # OS
            z = [zpair[0],zpair[1]] if zpair[0].pt()>zpair[1].pt() else [zpair[1],zpair[0]]
            if not bestZ: bestZ = z
            if z[0].pt()>bestZ[0].pt():
                bestZ = z
            if z[0].pt()==bestZ[0].pt() and z[1].pt()>bestZ[1].pt():
                bestZ = z
            #z = DiCandidate(*zpair[:2])
            #massdiff = abs(z.M()-ZMASS)
            #if massdiff<bestMassdiff:
            #    bestZ = zpair
            #    bestMassdiff = massdiff

        if not bestZ: return candidate # need a z candidate
        zcand = DiCandidate(*bestZ)
        if zcand.M()<60 or zcand.M()>120: return candidate

        tleps = [t for t in tleps if deltaR(t.eta(),t.phi(),zcand[0].eta(),zcand[0].phi())>0.8 and deltaR(t.eta(),t.phi(),zcand[1].eta(),zcand[1].phi())>0.8]
        if len(tleps)<1: return candidate

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['t'] = tleps[0]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['w'] = MetCompositeCandidate(self.pfmet,tleps[0])

        medLeps = self.getPassingCands('Medium',self.electrons,self.muons,self.taus)
        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps,0.4)

        # match jet to obj
        for c in ['t']:
            dr = 999
            j = None
            for jet in self.jets:
                jt = DiCandidate(jet,candidate[c])
                if jt.deltaR()<dr:
                    j = jet
                    dr = jt.deltaR()
            if j:
                candidate['{}jet'.format(c)] = j

        return candidate

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2','t']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,cands):
        nz = len(self.getPassingCands('Medium',self.electrons,self.muons))
        nt = len(self.getPassingCands('LooseNew',self.taus))
        return nz==2 and nt>0

    def trigger(self,cands):
        # accept MC, check trigger for data
        isData = self.event.isData()>0.5
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
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'SingleMuon'     : [
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

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('data'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='zTauFakeRateTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    zTauFakeRateAnalysis = ZTauFakeRateAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ZTauFakeRateTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       zTauFakeRateAnalysis.analyze()
       zTauFakeRateAnalysis.finish()
    except KeyboardInterrupt:
       zTauFakeRateAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
