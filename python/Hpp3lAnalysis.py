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

logger = logging.getLogger("Hpp3lAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class Hpp3lAnalysis(AnalysisBase):
    '''
    Hpp3l analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','hpp3lTree.root')
        outputTreeName = kwargs.pop('outputTreeName','Hpp3lTree')
        self.preselection = 'electrons_count+muons_count+taus_count>2'
        super(Hpp3lAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.threeLoose,'metFilter')
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.vetoFourth,'noFourthTightLepton')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',4])
        self.tree.add(self.getGenChannelString, 'genChannel', ['C',5])
        self.tree.add(self.getWZChannelString, 'wzChannel', ['C',4])

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
            self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.DoubleMediumIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg', 'I')
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
            #self.tree.add(lambda cands: self.event.Mu45_eta2p1Pass(), 'pass_Mu45_eta2p1', 'I')
            #self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg', 'I')

        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[0], 'triggerEfficiencyMC', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[1], 'triggerEfficiencyMCUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[2], 'triggerEfficiencyMCDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[0], 'triggerEfficiencyData', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[1], 'triggerEfficiencyDataUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[2], 'triggerEfficiencyDataDown', 'F')

        # 3 lepton
        self.addComposite('3l')
        self.addCompositeMet('3lmet')

        # hpp leptons
        self.addDiLepton('hpp')
        self.addCompositeMet('hppmet')
        self.addLepton('hpp1',doId=True,doScales=True,doFakes=True,doErrors=True)
        self.addLepton('hpp2',doId=True,doScales=True,doFakes=True,doErrors=True)

        # hm lepton
        self.addLeptonMet('hm')
        self.addLepton('hm1',doId=True,doScales=True,doFakes=True,doErrors=True)

        # wrong combination
        self.addDiLepton('hm1_hpp1')
        self.addDiLepton('hm1_hpp2')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('z1')
        self.addLepton('z2')

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('w1')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select 3l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'hpp1' : None,
            'hpp2' : None,
            'hm1' : None,
            'hpp' : None,
            'hppmet' : None,
            'hm': None,
            '3l': None,
            '3lmet': None,
            'hm1_hpp1': None,
            'hm1_hpp2': None,
            'z1' : Candidate(None),
            'z2' : Candidate(None),
            'w1' : Candidate(None),
            'z' : DiCandidate(Candidate(None),Candidate(None)),
            'w' : MetCompositeCandidate(self.pfmet,Candidate(None)),
            #'leadJet' : (),
            #'subleadJet' : (),
            'met': self.pfmet,
            'cleanJets' : [],
        }

        # get leptons
        leps = self.getPassingCands('Loose',self.electrons,self.muons,self.taus)
        medLeps = self.getPassingCands('Medium',self.electrons,self.muons,self.taus)
        if len(leps)<3: return candidate # need at least 3 leptons
        if len(medLeps)>3: return candidate # cant have more than 3 medium leptons

        # require ++- or --+
        if abs(sum([c.charge() for c in leps]))!=1: return candidate

        # get the candidates
        hppCands = []
        for trio in itertools.permutations(leps,3):
            # require ++-/--+
            if trio[0].charge()!=trio[1].charge(): continue
            if trio[0].charge()==trio[2].charge(): continue
            # require deltaR seperation of 0.02
            keep = True
            for i,j in itertools.combinations(range(3),2):
                dicand = DiCandidate(trio[i],trio[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            if not keep: continue
            # require lead e/m pt > 30/20, if all taus, lead 2 pt>40
            ems = [cand for cand in trio if cand.collName in ['electrons','muons']]
            ts = [cand for cand in trio if cand.collName in ['taus']]
            if len(ems)>0:
                pts_ems = sorted([cand.pt() for cand in ems])
                if pts_ems[-1]<30.: continue
                if len(ems)>1:
                    if pts_ems[-2]<20.: continue
            else:
                pts_ts = sorted([cand.pt() for cand in ts])
                if pts_ts[-2]<40.: continue
            ## allow at most 1 fake tau (Z)
            #numFake = 0
            #for t in ts:
            #    if t in leps and t not in medLeps: numFake += 1
            #if numFake>1: continue
            # its a good candidate
            hppCands += [trio]
        if not hppCands: return candidate
        hppCand = hppCands[0][:2]
        hmCand = hppCands[0][2]

        hpp1 = hppCand[0] if hppCand[0].pt()>hppCand[1].pt() else hppCand[1]
        hpp2 = hppCand[1] if hppCand[0].pt()>hppCand[1].pt() else hppCand[0]
        hm1 = hmCand

        candidate['hpp1'] = hpp1
        candidate['hpp2'] = hpp2
        candidate['hm1'] = hm1
        candidate['hpp'] = DiCandidate(hpp1,hpp2)
        candidate['hm'] = MetCompositeCandidate(self.pfmet,hm1)
        candidate['3l'] = CompositeCandidate(hpp1,hpp2,hm1)
        candidate['3lmet'] = MetCompositeCandidate(self.pfmet,hpp1,hpp2,hm1)
        candidate['hppmet'] = MetCompositeCandidate(self.pfmet,hpp1,hpp2)
        candidate['hm1_hpp1'] = DiCandidate(hm1,hpp1)
        candidate['hm1_hpp2'] = DiCandidate(hm1,hpp2)

        # add jet
        #jets = self.getCands(self.jets, lambda cand: cand.isLoose()>0.5)
        #if len(jets)==1:
        #    candidate['leadJet'] = jets[0]
        #    candidate['subleadJet'] = ('jets',-1)
        #if len(jets)>1:
        #    candidate['leadJet'] = jets[0]
        #    candidate['subleadJet'] = jets[1]
        #else:
        #    candidate['leadJet'] = ('jets',-1)
        #    candidate['subleadJet'] = ('jets',-1)

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for zpair in itertools.combinations(leps,3):
            if zpair[0].collName!=zpair[1].collName: continue # SF
            if zpair[0].charge()==zpair[1].charge(): continue # OS
            z = DiCandidate(*zpair[:2])
            massdiff = abs(z.M()-ZMASS)
            if massdiff<bestMassdiff:
                bestZ = zpair
                bestMassdiff = massdiff

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps+[hpp1,hpp2,hm1],0.4)

        if not bestZ: return candidate # no z

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        w = bestZ[2]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['w1'] = w
        candidate['w'] = MetCompositeCandidate(self.pfmet,w)

        return candidate

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['hpp1','hpp2','hm1']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    def getWZChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2','w1']:
            if not cands[c]: return 'aaa'
            chanString += self.getCollectionString(cands[c])
        return chanString

    def getGenChannelString(self,cands):
        '''Get the gen h++ channel'''
        chanString = ''
        pdgMap = {
            11: 'e',
            13: 'm',
            15: 't',
        }
        if 'HPlusPlusHMinusMinusHTo4L' in self.fileNames[0]: # h++h-- signal sample
            for s in [-1,1]:
                h = -1*s*9900041                     # h++ in pythia8
                for l1 in [s*11, s*13, s*15]:        # lepton 1
                    for l2 in [s*11, s*13, s*15]:    # lepton 2
                        if abs(l2)<abs(l1): continue # skip double counting
                        hasDecay = self.findDecay(h,l1,l2)
                        if hasDecay:
                            chanString += pdgMap[abs(l1)]
                            chanString += pdgMap[abs(l2)]
        elif 'HPlusPlusHMinusMinusHRTo4L' in self.fileNames[0]: # hr++hr-- signal sample
            for s in [-1,1]:
                h = -1*s*9900042                     # h++ in pythia8
                for l1 in [s*11, s*13, s*15]:        # lepton 1
                    for l2 in [s*11, s*13, s*15]:    # lepton 2
                        if abs(l2)<abs(l1): continue # skip double counting
                        hasDecay = self.findDecay(h,l1,l2)
                        if hasDecay:
                            chanString += pdgMap[abs(l1)]
                            chanString += pdgMap[abs(l2)]
        elif 'HPlusPlusHMinusHTo3L' in self.fileNames[0]: # h++h- signal sample
            for s in [-1,1]:
                h = -1*s*9900041                     # h++ in pythia8
                for l1 in [s*11, s*13, s*15]:        # lepton 1
                    for l2 in [s*11, s*13, s*15]:    # lepton 2
                        if abs(l2)<abs(l1): continue # skip double counting
                        hasDecay = self.findDecay(h,l1,l2)
                        if hasDecay:
                            chanString += pdgMap[abs(l1)]
                            chanString += pdgMap[abs(l2)]
            for s in [-1,1]:
                h = s*37                             # h+ in pythia8
                for l1 in [-1*s*11, -1*s*13, -1*s*15]:        # lepton 1
                    for l2 in [s*12, s*14, s*16]:    # neutrino
                        hasDecay = self.findDecay(h,l1,l2)
                        if hasDecay:
                            chanString += pdgMap[abs(l1)]
        else:
            chanString = 'a'
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,cands):
        return len(self.getPassingCands('Loose',self.electrons,self.muons,self.taus))>=3

    def vetoFourth(self,cands):
        return len(self.getPassingCands('Medium',self.electrons,self.muons,self.taus))<=3

    def trigger(self,cands):
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
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL',
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
                'Tau' : [
                    'DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg',
                ],
            }
            if isData and self.event.run()>=281639:
                triggerNames['MuonEG'] = [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
                    'Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ',
                ]

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
            'Tau',
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['hpp1','hpp2','hm1']]
        numTaus = [c.collName for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Ele27Tight','Mu17Mu8','Ele23Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','DoublePFTau35'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Ele27Tight','DoublePFTau35']
        else:
            triggerList = ['DoublePFTau35']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)




def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('hpp3l'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='hpp3lTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    hpp3lAnalysis = Hpp3lAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='Hpp3lTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       hpp3lAnalysis.analyze()
       hpp3lAnalysis.finish()
    except KeyboardInterrupt:
       hpp3lAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
