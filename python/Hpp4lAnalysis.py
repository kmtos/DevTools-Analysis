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

logger = logging.getLogger("Hpp4lAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class Hpp4lAnalysis(AnalysisBase):
    '''
    Hpp4l analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','hpp4lTree.root')
        outputTreeName = kwargs.pop('outputTreeName','Hpp4lTree')
        self.preselection = 'muons_count+electrons_count+taus_count>3'
        super(Hpp4lAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.fourLoose,'fourLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',5])
        self.tree.add(self.getGenChannelString, 'genChannel', ['C',5])
        self.tree.add(self.getZChannelString, 'zChannel', ['C',3])

        # event counts
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

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

        ## vbf
        #self.addJet('leadJet')
        #self.addJet('subleadJet')
        #self.addDiJet('dijet','leadJet','subleadJet')
        #self.tree.add(lambda cands: self.numCentralJets(cands,'isLoose',30), 'dijet_numCentralJetsLoose30', 'I')
        #self.tree.add(lambda cands: self.numCentralJets(cands,'isTight',30), 'dijet_numCentralJetsTight30', 'I')

        # 4 lepton
        self.addComposite('4l')
        self.addCompositeMet('4lmet')

        # hpp leptons
        self.addDiLepton('hpp')
        self.addCompositeMet('hppmet')
        self.addLepton('hpp1')
        self.tree.add(lambda cands: self.passMedium(cands['hpp1']),        'hpp1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hpp1']),         'hpp1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hpp1'])[0],     'hpp1_looseScale', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hpp1'])[1],     'hpp1_looseScaleUp', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hpp1'])[2],     'hpp1_looseScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp1'])[0],    'hpp1_mediumScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp1'])[1],    'hpp1_mediumScaleUp', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp1'])[2],    'hpp1_mediumScaleDown', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp1'])[0],     'hpp1_tightScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp1'])[1],     'hpp1_tightScaleUp', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp1'])[2],     'hpp1_tightScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp1'])[0], 'hpp1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp1'])[1], 'hpp1_mediumFakeRateUp', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp1'])[2], 'hpp1_mediumFakeRateDown', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp1'])[0],  'hpp1_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp1'])[1],  'hpp1_tightFakeRateUp', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp1'])[2],  'hpp1_tightFakeRateDown', 'F')
        self.addLepton('hpp2')
        self.tree.add(lambda cands: self.passMedium(cands['hpp2']),        'hpp2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hpp2']),         'hpp2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hpp2'])[0],     'hpp2_looseScale', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hpp2'])[1],     'hpp2_looseScaleUp', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hpp2'])[2],     'hpp2_looseScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp2'])[0],    'hpp2_mediumScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp2'])[1],    'hpp2_mediumScaleUp', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp2'])[2],    'hpp2_mediumScaleDown', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp2'])[0],     'hpp2_tightScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp2'])[1],     'hpp2_tightScaleUp', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp2'])[2],     'hpp2_tightScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp2'])[0], 'hpp2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp2'])[1], 'hpp2_mediumFakeRateUp', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp2'])[2], 'hpp2_mediumFakeRateDown', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp2'])[0],  'hpp2_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp2'])[1],  'hpp2_tightFakeRateUp', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp2'])[2],  'hpp2_tightFakeRateDown', 'F')

        # hmm leptons
        self.addDiLepton('hmm')
        self.addCompositeMet('hmmmet')
        self.addLepton('hmm1')
        self.tree.add(lambda cands: self.passMedium(cands['hmm1']),        'hmm1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hmm1']),         'hmm1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hmm1'])[0],     'hmm1_looseScale', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hmm1'])[1],     'hmm1_looseScaleUp', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hmm1'])[2],     'hmm1_looseScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm1'])[0],    'hmm1_mediumScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm1'])[1],    'hmm1_mediumScaleUp', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm1'])[2],    'hmm1_mediumScaleDown', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm1'])[0],     'hmm1_tightScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm1'])[1],     'hmm1_tightScaleUp', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm1'])[2],     'hmm1_tightScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm1'])[0], 'hmm1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm1'])[1], 'hmm1_mediumFakeRateUp', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm1'])[2], 'hmm1_mediumFakeRateDown', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm1'])[0],  'hmm1_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm1'])[1],  'hmm1_tightFakeRateUp', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm1'])[2],  'hmm1_tightFakeRateDown', 'F')
        self.addLepton('hmm2')
        self.tree.add(lambda cands: self.passMedium(cands['hmm2']),        'hmm2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hmm2']),         'hmm2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hmm2'])[0],     'hmm2_looseScale', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hmm2'])[1],     'hmm2_looseScaleUp', 'F')
        self.tree.add(lambda cands: self.looseScale(cands['hmm2'])[2],     'hmm2_looseScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm2'])[0],    'hmm2_mediumScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm2'])[1],    'hmm2_mediumScaleUp', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm2'])[2],    'hmm2_mediumScaleDown', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm2'])[0],     'hmm2_tightScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm2'])[1],     'hmm2_tightScaleUp', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm2'])[2],     'hmm2_tightScaleDown', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm2'])[0], 'hmm2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm2'])[1], 'hmm2_mediumFakeRateUp', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm2'])[2], 'hmm2_mediumFakeRateDown', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm2'])[0],  'hmm2_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm2'])[1],  'hmm2_tightFakeRateUp', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm2'])[2],  'hmm2_tightFakeRateDown', 'F')

        # wrong combination
        self.addDiLepton('hmm1_hpp1')
        self.addDiLepton('hmm1_hpp2')
        self.addDiLepton('hmm2_hpp1')
        self.addDiLepton('hmm2_hpp2')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('z1')
        self.addLepton('z2')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'hpp1' : None,
            'hpp2' : None,
            'hmm1' : None,
            'hmm2' : None,
            'hpp' : None,
            'hppmet' : None,
            'hmm' : None,
            'hmmmet' : None,
            '4l' : None,
            '4lmet' : None,
            'hmm1_hpp1' : None,
            'hmm1_hpp2' : None,
            'hmm2_hpp1' : None,
            'hmm2_hpp2' : None,
            'z1' : Candidate(None),
            'z2' : Candidate(None),
            'z' : DiCandidate(Candidate(None),Candidate(None)),
            #'leadJet' : (),
            #'subleadJet' : (),
            'met': self.pfmet,
            'cleanJets' : [],
        }

        # get leptons
        leps = self.getPassingCands('Loose',self.electrons,self.muons,self.taus)
        medLeps = self.getPassingCands('Medium',self.electrons,self.muons,self.taus)
        if len(leps)<4: return candidate # need at least 4 leptons


        # get the candidates
        hppCands = []
        for quad in itertools.permutations(leps,4):
            # require ++--
            if quad[0].charge()+quad[1].charge()!=2: continue
            if quad[2].charge()+quad[3].charge()!=-2: continue
            # require deltaR seperation of 0.02 m(ll)>12
            keep = True
            for i,j in itertools.combinations(range(4),2):
                dicand = DiCandidate(quad[i],quad[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            if not keep: continue
            # require lead e/m pt > 25, if all taus, lead 2 pt>40
            ems = [cand for cand in quad if cand.collName in ['electrons','muons']]
            ts = [cand for cand in quad if cand.collName in ['taus']]
            if len(ems)>0:
                pts_ems = sorted([cand.pt() for cand in ems])
                if pts_ems[-1]<30.: continue
                if len(ems)>1:
                    if pts_ems[-2]<20.: continue
            else:
                pts_ts = [cand.pt() for cand in ts]
                if sorted(pts_ts)[-2]<40.: continue
            ## allow at most 2 fake taus (Z)
            #numFake = 0
            #for t in ts:
            #    if t in leps and t not in medLeps: numFake += 1
            #if numFake>2: continue
            # its a good candidate
            hppCands += [quad]
        if not hppCands: return candidate

        # sort by highest st
        highestSt = 0
        bestCand = []
        for quad in hppCands:
            st = sum([cand.pt() for cand in quad])
            if st>highestSt:
                bestCand = quad
                highestSt = st

        hpp1 = bestCand[0] if bestCand[0].pt()>bestCand[1].pt() else bestCand[1]
        hpp2 = bestCand[1] if bestCand[0].pt()>bestCand[1].pt() else bestCand[0]
        hmm1 = bestCand[2] if bestCand[2].pt()>bestCand[3].pt() else bestCand[3]
        hmm2 = bestCand[3] if bestCand[2].pt()>bestCand[3].pt() else bestCand[2]

        candidate['hpp1'] = hpp1
        candidate['hpp2'] = hpp2
        candidate['hmm1'] = hmm1
        candidate['hmm2'] = hmm2
        candidate['hpp'] = DiCandidate(hpp1,hpp2)
        candidate['hppmet'] = MetCompositeCandidate(self.pfmet,hpp1,hpp2)
        candidate['hmm'] = DiCandidate(hmm1,hmm2)
        candidate['hmmmet'] = MetCompositeCandidate(self.pfmet,hmm1,hmm2)
        candidate['4l'] = CompositeCandidate(hpp1,hpp2,hmm1,hmm2)
        candidate['4lmet'] = MetCompositeCandidate(self.pfmet,hpp1,hpp2,hmm1,hmm2)
        candidate['hmm1_hpp1'] = DiCandidate(hmm1,hpp1)
        candidate['hmm1_hpp2'] = DiCandidate(hmm1,hpp2)
        candidate['hmm2_hpp1'] = DiCandidate(hmm2,hpp1)
        candidate['hmm2_hpp2'] = DiCandidate(hmm2,hpp2)

        # get invariant masses
        bestZ = ()
        bestMassdiff = 99999
        for zzquad in itertools.combinations(leps,4):
            if zzquad[0].collName!=zzquad[1].collName: continue # SF
            if zzquad[0].charge()==zzquad[1].charge(): continue # OS
            z = DiCandidate(*zzquad[:2])
            massdiff = abs(z.M()-ZMASS)
            if massdiff<bestMassdiff:
                bestZ = zzquad
                bestMassdiff = massdiff

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps+[hpp1,hpp2,hmm1,hmm2],0.4)

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = [bestZ[0],bestZ[1]] if bestZ[0].pt()>bestZ[1].pt() else [bestZ[1],bestZ[0]]
        candidate['z1'] = z[0]
        candidate['z2'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])

        return candidate


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['hpp1','hpp2','hmm1','hmm2']:
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

    def getZChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2']:
            if not cands[c]: return 'aa'
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def fourLoose(self,cands):
        return len(self.getPassingCands('Loose',self.electrons,self.muons,self.taus))>=4

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
        candList = [cands[c] for c in ['hpp1','hpp2','hmm1','hmm2']]
        numTaus = [c.collName for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Ele27Tight','Mu17Mu8','Ele23Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12','DoublePFTau35'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Ele27Tight','Mu17Mu8','Ele23Ele12','DoublePFTau35']
        elif numTaus==3:
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

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('hpp'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='hpp4lTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    hpp4lAnalysis = Hpp4lAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='Hpp4lTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       hpp4lAnalysis.analyze()
       hpp4lAnalysis.finish()
    except KeyboardInterrupt:
       hpp4lAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
