# Hpp3lAnalysis.py
# for hpp3l analysis

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight

from Candidates import *

import sys
import itertools
import operator

import ROOT

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
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.vetoFourth,'noFourthTightLepton')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',4])
        self.tree.add(self.getGenChannelString, 'genChannel', ['C',5])
        self.tree.add(self.getWZChannelString, 'wzChannel', ['C',4])

        # event counts
        self.tree.add(lambda cands: self.numJets('isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets('isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets('passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
        self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        ## vbf
        #self.addJet('leadJet')
        #self.addJet('subleadJet')
        #self.addDiJet('dijet')
        #self.tree.add(lambda cands: self.numCentralJets(rcands,'isLoose',30), 'dijet_numCentralJetsLoose30', 'I')
        #self.tree.add(lambda cands: self.numCentralJets(rcands,'isTight',30), 'dijet_numCentralJetsTight30', 'I')

        # 3 lepton
        self.addComposite('3l')
        self.addCompositeMet('3lmet')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hpp1'],cands['hpp2'],cands['hm1']), '3l_zeppenfeld','F')

        # hpp leptons
        self.addDiLepton('hpp')
        self.addCompositeMet('hppmet')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hpp1'],cands['hpp2']), 'hpp_zeppenfeld','F')
        self.addLepton('hpp1')
        self.tree.add(lambda cands: self.passMedium(cands['hpp1']), 'hpp1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hpp1']), 'hpp1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hpp1']), 'hpp1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp1']), 'hpp1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp1']), 'hpp1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp1']), 'hpp1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp1']), 'hpp1_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hpp1']), 'hpp1_zeppenfeld','F')
        self.addLepton('hpp2')
        self.tree.add(lambda cands: self.passMedium(cands['hpp2']), 'hpp2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hpp2']), 'hpp2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hpp2']), 'hpp2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hpp2']), 'hpp2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hpp2']), 'hpp2_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hpp2']), 'hpp2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hpp2']), 'hpp2_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hpp2']), 'hpp2_zeppenfeld','F')

        # hm lepton
        self.addLeptonMet('hm')
        self.addLepton('hm1')
        self.tree.add(lambda cands: self.passMedium(cands['hm1']), 'hm1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hm1']), 'hm1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hm1']), 'hm1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hm1']), 'hm1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hm1']), 'hm1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hm1']), 'hm1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hm1']), 'hm1_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hm1']), 'hm1_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('hm1_hpp1')
        self.addDiLepton('hm1_hpp2')

        # z leptons
        self.addDiLepton('z')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z1'],cands['z2']), 'z_zeppenfeld','F')
        self.addLepton('z1')
        self.tree.add(lambda cands: self.passMedium(cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z1']), 'z1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z1']), 'z1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z1']), 'z1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z1']), 'z1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z1']), 'z1_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z1']), 'z1_zeppenfeld','F')
        self.addLepton('z2')
        self.tree.add(lambda cands: self.passMedium(cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z2']), 'z2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z2']), 'z2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z2']), 'z2_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z2']), 'z2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z2']), 'z2_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z2']), 'z2_zeppenfeld','F')

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('w1')
        self.tree.add(lambda cands: self.passMedium(cands['w1']), 'w1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['w1']), 'w1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['w1']), 'w1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['w1']), 'w1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['w1']), 'w1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['w1']), 'w1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['w1']), 'w1_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['w1']), 'w1_zeppenfeld','F')

        # met
        self.addMet('met')

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
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        medLeps = self.getPassingCands('Medium')
        if len(leps)<3: return candidate # need at least 3 leptons
        if len(medLeps)>3: return candidate # cant have more than 3 medium leptons

        # require ++- or --+
        if abs(sum([c.charge() for c in leps]))!=1: return candidate

        # get the candidates
        hppCands = []
        for trio in itertools.permutations(leps,3):
            # require ++-/--+
            if trio[0].charge()!=trio[1].charge(): continue
            # require deltaR seperation of 0.02
            keep = True
            for i,j in itertools.combinations(range(3),2):
                dicand = DiCandidate(trio[i],trio[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            if not keep: continue
            # require lead e/m pt > 20, if all taus, lead 2 pt>35
            ems = [cand for cand in trio if cand.collName in ['electrons','muons']]
            ts = [cand for cand in trio if cand.collName in ['taus']]
            if len(ems)>0:
                pts_ems = [cand.pt() for cand in ems]
                if max(pts_ems)<20.: continue
            else:
                pts_ts = [cand.pt() for cand in ts]
                if sorted(pts_ts)[-2]<35.: continue
            # allow at most 1 fake tau (Z)
            numFake = 0
            for t in ts:
                if t in leps and t not in medLeps: numFake += 1
            if numFake>1: continue
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
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDLooseIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedVeto',cand)
        else:
            return 1.

    def mediumScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedMedium',cand)
        else:
            return 1.

    def tightScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedTight',cand)
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
        for coll in [self.electrons,self.muons,self.taus]:
            cands += self.getCands(coll,passMode)
        return cands

    def numJets(self,mode,pt):
        return len(
            self.getCands(
                self.jets,
                lambda cand: getattr(cand,mode)()>0.5 and cand.pt()>pt
            )
        )

    def numCentralJets(self,cands,mode,pt):
        if not cands['leadJet']: return -1
        if not cands['subleadJet']: return -1
        eta1 = cands['leadJet'].eta()
        eta2 = cands['subleadJet'].eta()
        mineta = min(eta1,eta2)
        maxeta = max(eta1,eta2)
        return len(
            self.getCands(
                self.jets,
                lambda cand: getattr(cand,mode)()>0.5
                             and cand.pt()>pt
                             and cand.eta()>mineta
                             and cand.eta()<maxeta
            )
        )

    def zeppenfeld(self,cands,*probeCands):
        if not cands['leadJet']: return -10.
        if not cands['subleadJet']: return -10.
        eta1 = cands['leadJet'].eta()
        eta2 = cands['subleadJet'].eta()
        meaneta = (eta1+eta2)/2
        composite = CompositeCandidate(*probeCands)
        eta = composite.eta()
        return eta-meaneta


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
        else:
            chanString = 'a'
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,cands):
        return len(self.getPassingCands('Loose'))>=3

    def vetoFourth(self,cands):
        return len(self.getPassingCands('Medium'))<=3

    def trigger(self,cands):
        # accept MC, check trigger for data
        if self.event.isData()<0.5: return True
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
        candList = [cands[c] for c in ['hpp1','hpp2','hm1']]
        numTaus = [c.collName for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','DoublePFTau35']
        else:
            triggerList = ['DoublePFTau35']
        return self.triggerScales.getDataEfficiency(triggerList,candList)







