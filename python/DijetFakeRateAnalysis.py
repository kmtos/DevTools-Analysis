from AnalysisBase import AnalysisBase
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from utilities import ZMASS, deltaPhi, deltaR
from Candidates import *

import sys
import itertools
import operator

sys.argv.append('-b')
import ROOT
sys.argv.pop()

class DijetFakeRateAnalysis(AnalysisBase):
    '''
    Select a single lepton to perform a dijet control fake rate
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','dijetFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','DijetFakeRateTree')
        super(DijetFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.vetoSecond,'vetoSecond')
        self.cutTree.add(self.metVeto,'metVeto')
        self.cutTree.add(self.mtVeto,'mtVeto')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',2])

        # event counts
        self.tree.add(lambda cands: self.numJets('isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets('isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets('passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passMedium)), 'numMediumElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passMedium)), 'numMediumMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVLPass(), 'pass_Mu8_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Mu24_TrkIsoVVLPass(), 'pass_Mu24_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Mu34_TrkIsoVVLPass(), 'pass_Mu34_TrkIsoVVL', 'I')
        self.tree.add(lambda cands: self.event.Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Ele23_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele23_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerPrescale, 'triggerPrescale', 'F')

        # lead jet
        self.addJet('leadJet')

        # lepton
        self.addLeptonMet('w')
        self.addLepton('l1')
        self.tree.add(lambda cands: self.passMedium(cands['l1']), 'l1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['l1']), 'l1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['l1']), 'l1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['l1']), 'l1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['l1']), 'l1_tightScale', 'F')

        # met
        self.addMet('met')

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self):
        candidate = {
            'l1' : None,
            'w' : None,
            'leadJet' : None,
            'met': self.pfmet
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        if len(leps)<1: return candidate # need at least 1 lepton

        # choose highest pt
        bestCand = 0
        bestPt = 0
        for l in leps:
           if l.pt()>bestPt:
               bestPt = l.pt()
               bestCand = l

        # add jet
        jets = self.getCands(self.jets, lambda cand: cand.isLoose()>0.5)
        if len(jets)<0: return candidate # need a recoil jet

        # choose highest pt
        bestJet = 0
        bestPt = 0
        for j in jets:
           if j.pt()>bestPt:
               bestPt = j.pt()
               bestJet = j

        candidate['l1'] = bestCand
        candidate['w'] = MetCompositeCandidate(self.pfmet,bestCand)
        candidate['leadJet'] = bestJet

        return candidate

    ##################
    ### lepton IDs ###
    ##################
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

    def numJets(self,mode,pt):
        return len(
            self.getCands(
                self.jets,
                lambda cand: getattr(cand,mode)()>0.5 and cand.pt()>pt
            )
        )

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['l1']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def vetoSecond(self,cands):
        return len(self.getPassingCands('Loose'))==1

    def metVeto(self,cands):
        return cands['met'].et()<25

    def mtVeto(self,cands):
        if not cands['w']: return False
        return cands['w'].Mt()<25

    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
        if not cands['l1']: return False
        triggerNames = {
            'DoubleMuon'     : [
                ['Mu8_TrkIsoVVL', 0],
                ['Mu17_TrkIsoVVL', 20],
                #['Mu24_TrkIsoVVL', 30],
                #['Mu34_TrkIsoVVL', 40],
            ],
            'DoubleEG'       : [
                ['Ele12_CaloIdL_TrackIdL_IsoVL', 0],
                ['Ele17_CaloIdL_TrackIdL_IsoVL', 20],
                #['Ele23_CaloIdL_TrackIdL_IsoVL', 30],
            ],
        }
        # here we need to accept only a specific trigger after a certain pt threshold
        pt = cands['l1'].pt()
        dataset = 'DoubleEG' if isinstance(cands['l1'],Electron) else 'DoubleMuon'
        # accept the event only if it is triggered in the current dataset
        reject = True if self.event.isData()>0.5 else False
        if dataset in self.fileNames[0]: reject = False
        # now pick the appropriate trigger for the pt
        theTrigger = ''
        for trig, ptThresh in triggerNames[dataset]:
            if pt < ptThresh: break
            theTrigger = trig
        # and check if we pass
        passTrigger = getattr(self.event,'{0}Pass'.format(theTrigger))()
        if passTrigger>0.5: return False if reject else True
        return False

    def triggerEfficiency(self,cands):
        candList = [cands['l1']]
        # select via pt and flavor
        pt = cands['l1'].pt()
        if isinstance(cands['l1'],Electron):
            if pt<20:
                triggerList = ['Ele17_Ele12Leg2']
            else:
                triggerList = ['Ele17_Ele12Leg1']
        else:
            if pt<20:
                triggerList = ['Mu17_Mu8Leg2']
            else:
                triggerList = ['Mu17_Mu8Leg1']
        return self.triggerScales.getDataEfficiency(triggerList,candList)

    def triggerPrescale(self,cands):
        # select via pt and flavor
        pt = cands['l1'].pt()
        if isinstance(cands['l1'],Electron):
            if pt<20:
                trigger = 'Ele17_Ele12Leg2'
            else:
                trigger = 'Ele17_Ele12Leg1'
        else:
            if pt<20:
                trigger = 'Mu17_Mu8Leg2'
            else:
                trigger = 'Mu17_Mu8Leg1'
        return self.triggerPrescales.getPrescale(trigger)

