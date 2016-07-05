# TauChargeAnalysis.py
# for DY analysis

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight

from Candidates import *

import itertools
import operator

import ROOT

class TauChargeAnalysis(AnalysisBase):
    '''
    TauCharge analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','tauChargeTree.root')
        outputTreeName = kwargs.pop('outputTreeName','TauChargeTree')
        super(TauChargeAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.twoLoose,'twoLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu22Pass(), 'pass_IsoMu22', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu22Pass(), 'pass_IsoTkMu22', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # z leptons
        self.addDiLepton('z')
        self.addLepton('m')
        self.tree.add(lambda cands: self.passMedium(cands['m']), 'm_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['m']), 'm_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['m']), 'm_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['m']), 'm_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['m']), 'm_tightScale', 'F')
        self.addLepton('t')
        self.tree.add(lambda cands: self.passMedium(cands['t']), 't_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['t']), 't_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['t']), 't_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['t']), 't_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['t']), 't_tightScale', 'F')

        # w lepton
        self.addLeptonMet('wm')
        self.addLeptonMet('wt')

        # met
        self.addMet('met')

    ############################
    ### select DY candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'm' : None,
            't' : None,
            'z' : None,
            'wm': None,
            'wt': None,
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        if len(leps)<2: return candidate # need at least 2 leptons

        # get invariant masses
        bestZ = ()
        bestst = 0
        for zpair in itertools.permutations(leps,2):
            if zpair[0].collName!='muons': continue
            if zpair[1].collName!='taus': continue
            if zpair[0].pt()<20: continue
            if zpair[1].pt()<20: continue
            z = DiCandidate(*zpair)
            if z.deltaR()<0.02: continue
            st = z.st()
            if st>bestst:
                bestZ = zpair
                bestst = st

        if not bestZ: return candidate # need a z candidate

        # and sort pt of Z
        z = bestZ
        candidate['m'] = z[0]
        candidate['t'] = z[1]
        candidate['z'] = DiCandidate(z[0],z[1])
        candidate['wm'] = MetCompositeCandidate(self.pfmet,z[0])
        candidate['wt'] = MetCompositeCandidate(self.pfmet,z[1])


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
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDLooseIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedVeto',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDVeto',cand)
        else:
            return 1.

    def mediumScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedMedium',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDMedium',cand)
        else:
            return 1.

    def tightScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedTight',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDTight',cand)
        else:
            return 1.

    def mediumFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppMedium','HppLoose')

    def tightFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppLoose')

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
        for coll in [self.muons,self.taus]:
            cands += self.getCands(coll,passMode)
        return cands


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['m','t']:
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
                'SingleMuon'     : [
                    'IsoMu20',
                    'IsoTkMu20',
                ],
            }
        else:
            triggerNames = {
                'SingleMuon'     : [
                    'IsoMu22',
                    'IsoTkMu22',
                ],
            }

        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
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
        candList = [cands['m']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu22ORIsoTkMu22']
        return self.triggerScales.getDataEfficiency(triggerList,candList)










