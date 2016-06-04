# Hpp4lAnalysis.py
# for hpp4l analysis

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight, passHZZLoose, passHZZTight

from Candidates import *

import sys
import itertools
import operator

import ROOT

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
        self.tree.add(lambda cands: self.numJets('isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numJets('isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda cands: self.numJets('passCSVv2T',30), 'numBjetsTight30', 'I')
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
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
        self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
        self.tree.add(lambda cands: self.event.DoubleMediumIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        ## vbf
        #self.addJet('leadJet')
        #self.addJet('subleadJet')
        #self.addDiJet('dijet','leadJet','subleadJet')
        #self.tree.add(lambda cands: self.numCentralJets(cands,'isLoose',30), 'dijet_numCentralJetsLoose30', 'I')
        #self.tree.add(lambda cands: self.numCentralJets(cands,'isTight',30), 'dijet_numCentralJetsTight30', 'I')

        # 4 lepton
        self.addComposite('4l')
        self.addCompositeMet('4lmet')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hpp1'],cands['hpp2'],cands['hmm1'],cands['hmm2']), '4l_zeppenfeld','F')

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

        # hmm leptons
        self.addDiLepton('hmm')
        self.addCompositeMet('hmmmet')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hmm1'],cands['hmm2']), 'hmm_zeppenfeld','F')
        self.addLepton('hmm1')
        self.tree.add(lambda cands: self.passMedium(cands['hmm1']), 'hmm1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hmm1']), 'hmm1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hmm1']), 'hmm1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm1']), 'hmm1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm1']), 'hmm1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm1']), 'hmm1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm1']), 'hmm1_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hmm1']), 'hmm1_zeppenfeld','F')
        self.addLepton('hmm2')
        self.tree.add(lambda cands: self.passMedium(cands['hmm2']), 'hmm2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['hmm2']), 'hmm2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['hmm2']), 'hmm2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['hmm2']), 'hmm2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['hmm2']), 'hmm2_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['hmm2']), 'hmm2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['hmm2']), 'hmm2_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['hmm2']), 'hmm2_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('hmm1_hpp1')
        self.addDiLepton('hmm1_hpp2')
        self.addDiLepton('hmm2_hpp1')
        self.addDiLepton('hmm2_hpp2')

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

        # met
        self.addMet('met')

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
            'z1' : Candidate(),
            'z2' : Candidate(),
            'z' : DiCandidate(Candidate(),Candidate()),
            #'leadJet' : (),
            #'subleadJet' : (),
            'met': self.pfmet,
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        medLeps = self.getPassingCands('Medium')
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
            # require lead e/m pt > 20, if all taus, lead 2 pt>35
            ems = [cand for cand in quad if cand.collName in ['electrons','muons']]
            ts = [cand for cand in quad if cand.collName in ['taus']]
            if len(ems)>0:
                pts_ems = [cand.pt() for cand in ems]
                if max(pts_ems)<20.: continue
            else:
                pts_ts = [cand.pt() for cand in ts]
                if sorted(pts_ts)[-2]<35.: continue
            # allow at most 2 fake taus (Z)
            numFake = 0
            for t in ts:
                if t in leps and t not in medLeps: numFake += 1
            if numFake>2: continue
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
        for zpair in itertools.combinations(leps,2):
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

        return candidate

    ##################
    ### lepton IDs ###
    ##################
    def passLoose(self,cand):
        return passHppLoose(cand)
        #return passHZZLoose(cand)

    def passMedium(self,cand):
        return passHppMedium(cand)
        #return passHZZTight(cand)

    def passTight(self,cand):
        return passHppTight(cand)
        #return passHZZTight(cand)

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
        return len(self.getPassingCands('Loose'))>=4

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
                    'Ele17_Ele12_CaloIdL_TrackIdL_IsoVL',
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
        candList = [cands[c] for c in ['hpp1','hpp2','hmm1','hmm2']]
        numTaus = [c.collName for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12','DoublePFTau35']
        elif numTaus==3:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','DoublePFTau35']
        else:
            triggerList = ['DoublePFTau35']
        return self.triggerScales.getDataEfficiency(triggerList,candList)









