from AnalysisBase import AnalysisBase
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight
from utilities import ZMASS, deltaPhi, deltaR

import sys
import itertools
import operator

sys.argv.append('-b')
import ROOT
sys.argv.pop()

class WTauFakeRateAnalysis(AnalysisBase):
    '''
    Select a muon + tau for fake rate in a W control region
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','WTauFakeRateTree.root')
        outputTreeName = kwargs.pop('outputTreeName','WTauFakeRateTree')
        super(WTauFakeRateAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',3])

        # event counts
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'taus',self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'taus',self.passMedium)),'numMediumTaus', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'taus',self.passTight)), 'numTightTaus', 'I')


        # trigger
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'IsoMu20Pass'), 'pass_IsoMu20', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'IsoTkMu20Pass'), 'pass_IsoTkMu20', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # lepton
        # mu tag
        self.addLeptonMet('wm','m',('pfmet',0))
        self.addLepton('m')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['m']), 'm_tightScale', 'F')

        # tau probe
        self.addLeptonMet('wt','t',('pfmet',0))
        self.addLepton('t')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['t']), 't_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['t']), 't_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['t']), 't_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['t']), 't_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['t']), 't_tightScale', 'F')

        # dilepton combination
        self.addDiLepton('z','m','t')

        # met
        self.addMet('met',('pfmet',0))

    #############################
    ### select fake candidate ###
    #############################
    def selectCandidates(self,rtrow):
        candidate = {
            'm' : (),
            't' : (),
        }

        # get leptons
        muons = self.getCands(rtrow,'muons',self.passTight)
        loosemuons = self.getCands(rtrow,'muons',self.passLoose)
        taus = self.getCands(rtrow,'taus',self.passLoose)
        if len(muons)!=1: return candidate # need 1 tight muon
        if len(loosemuons)>1: return candidate # dont allow another loose muon
        if len(taus)<1: return candidate # need at least 1 tau 

        mupt = self.getObjectVariable(rtrow,muons[0],'pt')
        if mupt<25: return candidate

        # select highest pt tau
        tauCand = None
        highest = 0
        for tau in taus:
            pt = self.getObjectVariable(rtrow,tau,'pt')
            if pt > highest:
                highest = pt
                tauCand = tau

        if not tauCand: return candidate

        taupt = self.getObjectVariable(rtrow,tauCand,'pt')
        if taupt<20: return candidate

        candidate['m'] = muons[0]
        candidate['t'] = tauCand

        return candidate

    ##################
    ### lepton IDs ###
    ##################
    def passLoose(self,rtrow,cand):
        return passHppLoose(self,rtrow,cand)

    def passMedium(self,rtrow,cand):
        return passHppMedium(self,rtrow,cand)

    def passTight(self,rtrow,cand):
        return passHppTight(self,rtrow,cand)

    def looseScale(self,rtrow,cand):
        if cand[0]=='muons':
            return self.leptonScales.getScale(rtrow,'MediumIDLooseIso',cand)
        elif cand[0]=='electrons':
            return self.leptonScales.getScale(rtrow,'CutbasedVeto',cand) # TODO, fix
        else:
            return 1.

    def mediumScale(self,rtrow,cand):
        if cand[0]=='muons':
            return self.leptonScales.getScale(rtrow,'MediumIDTightIso',cand)
        elif cand[0]=='electrons':
            return self.leptonScales.getScale(rtrow,'CutbasedMedium',cand)
        else:
            return 1.

    def tightScale(self,rtrow,cand):
        if cand[0]=='muons':
            return self.leptonScales.getScale(rtrow,'MediumIDTightIso',cand)
        elif cand[0]=='electrons':
            return self.leptonScales.getScale(rtrow,'CutbasedTight',cand)
        else:
            return 1.

    def getPassingCands(self,rtrow,mode):
        if mode=='Loose':
            passMode = self.passLoose
        elif mode=='Medium':
            passMode = self.passMedium
        elif mode=='Tight':
            passMode = self.passTight
        else:
            return []
        cands = []
        for coll in ['muons']:
            cands += self.getCands(rtrow,coll,passMode)
        return cands

    def numJets(self,rtrow,mode,pt):
        return len(
            self.getCands(
                rtrow,
                'jets',
                lambda rtrow,cand: self.getObjectVariable(rtrow,cand,mode)>0.5 
                                   and self.getObjectVariable(rtrow,cand,'pt')>pt
            )
        )

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,rtrow,cands):
        '''Get the channel string'''
        chanString = 'mt'
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,rtrow,cands):
        # accept MC, check trigger for data
        if rtrow.isData<0.5: return True
        triggerNames = {
            'SingleMuon'     : [
                'IsoMu20',
                'IsoTkMu20',
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
        reject = True if rtrow.isData>0.5 else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0]: reject = False
            for trigger in triggerNames[dataset]:
                var = '{0}Pass'.format(trigger)
                passTrigger = self.getTreeVariable(rtrow,var)
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0]: break
        return False

    def triggerEfficiency(self,rtrow,cands):
        candList = [cands['m']]
        triggerList = ['IsoMu20_OR_IsoTkMu20']
        return self.triggerScales.getDataEfficiency(rtrow,triggerList,candList)

