# Hpp3lAnalysis.py
# for hpp3l analysis

from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight

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
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'isLoose',30), 'numJetsLoose30', 'I')
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'isTight',30), 'numJetsTight30', 'I')
        self.tree.add(lambda rtrow,cands: self.numJets(rtrow,'passCSVv2T',30), 'numBjetsTight30', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'electrons',self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'electrons',self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'muons',self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda rtrow,cands: len(self.getCands(rtrow,'muons',self.passTight)), 'numTightMuons', 'I')

        # trigger
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass'), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass'), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass'), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass'), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass'), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'IsoMu20Pass'), 'pass_IsoMu20', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'IsoTkMu20Pass'), 'pass_IsoTkMu20', 'I')
        self.tree.add(lambda rtrow,cands: self.getTreeVariable(rtrow,'Ele23_WPLoose_GsfPass'), 'pass_Ele23_WPLoose_Gsf', 'I')
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')

        # vbf
        self.addJet('leadJet')
        self.addJet('subleadJet')
        self.addDiJet('dijet','leadJet','subleadJet')
        self.tree.add(lambda rtrow,cands: self.numCentralJets(rtrow,cands,'isLoose',30), 'dijet_numCentralJetsLoose30', 'I')
        self.tree.add(lambda rtrow,cands: self.numCentralJets(rtrow,cands,'isTight',30), 'dijet_numCentralJetsTight30', 'I')

        # 3 lepton
        self.addComposite('3l','hpp1','hpp2','hm1')
        self.addCompositeMet('3lmet',('pfmet',0),'hpp1','hpp2','hm1')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['hpp1'],cands['hpp2'],cands['hm1']), '3l_zeppenfeld','F')

        # hpp leptons
        self.addDiLepton('hpp','hpp1','hpp2')
        self.addCompositeMet('hppmet',('pfmet',0),'hpp1','hpp2')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['hpp1'],cands['hpp2']), 'hpp_zeppenfeld','F')
        self.addLepton('hpp1')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['hpp1']), 'hpp1_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['hpp1']), 'hpp1_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['hpp1']), 'hpp1_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['hpp1']), 'hpp1_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['hpp1']), 'hpp1_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['hpp1']), 'hpp1_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['hpp1']), 'hpp1_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['hpp1']), 'hpp1_zeppenfeld','F')
        self.addLepton('hpp2')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['hpp2']), 'hpp2_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['hpp2']), 'hpp2_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['hpp2']), 'hpp2_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['hpp2']), 'hpp2_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['hpp2']), 'hpp2_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['hpp2']), 'hpp2_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['hpp2']), 'hpp2_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['hpp2']), 'hpp2_zeppenfeld','F')

        # hm lepton
        self.addLeptonMet('hm','hm1',('pfmet',0))
        self.addLepton('hm1')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['hm1']), 'hm1_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['hm1']), 'hm1_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['hm1']), 'hm1_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['hm1']), 'hm1_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['hm1']), 'hm1_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['hm1']), 'hm1_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['hm1']), 'hm1_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['hm1']), 'hm1_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('hm1_hpp1','hm1','hpp1')
        self.addDiLepton('hm1_hpp2','hm1','hpp2')

        # z leptons
        self.addDiLepton('z','z1','z2')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['z1'],cands['z2']), 'z_zeppenfeld','F')
        self.addLepton('z1')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['z1']), 'z1_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['z1']), 'z1_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['z1']), 'z1_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['z1']), 'z1_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['z1']), 'z1_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['z1']), 'z1_zeppenfeld','F')
        self.addLepton('z2')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['z2']), 'z2_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['z2']), 'z2_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['z2']), 'z2_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['z2']), 'z2_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['z2']), 'z2_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['z2']), 'z2_zeppenfeld','F')

        # w lepton
        self.addLeptonMet('w','w1',('pfmet',0))
        self.addLepton('w1')
        self.tree.add(lambda rtrow,cands: self.passMedium(rtrow,cands['w1']), 'w1_passMedium', 'I')
        self.tree.add(lambda rtrow,cands: self.passTight(rtrow,cands['w1']), 'w1_passTight', 'I')
        self.tree.add(lambda rtrow,cands: self.looseScale(rtrow,cands['w1']), 'w1_looseScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumScale(rtrow,cands['w1']), 'w1_mediumScale', 'F')
        self.tree.add(lambda rtrow,cands: self.tightScale(rtrow,cands['w1']), 'w1_tightScale', 'F')
        self.tree.add(lambda rtrow,cands: self.mediumFakeRate(rtrow,cands['w1']), 'w1_mediumFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.tightFakeRate(rtrow,cands['w1']), 'w1_tightFakeRate', 'F')
        self.tree.add(lambda rtrow,cands: self.zeppenfeld(rtrow,cands,cands['w1']), 'w1_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('w1_z1','w1','z1')
        self.addDiLepton('w1_z2','w1','z2')

        # met
        self.addMet('met',('pfmet',0))

    ############################
    ### select 3l candidates ###
    ############################
    def selectCandidates(self,rtrow):
        candidate = {
            'hpp1' : (),
            'hpp2' : (),
            'hm1' : (),
            'z1' : (),
            'z2' : (),
            'w1' : (),
            'leadJet' : (),
            'subleadJet' : (),
        }

        # get leptons
        colls = ['electrons','muons','taus']
        pts = {}
        etas = {}
        phis = {}
        p4s = {}
        charges = {}
        leps = []
        leps = self.getPassingCands(rtrow,'Loose')
        medLeps = self.getPassingCands(rtrow,'Medium')
        if len(leps)<3: return candidate # need at least 3 leptons
        if len(medLeps)>3: return candidate # cant have more than 3 medium leptons

        for cand in leps:
            pts[cand] = self.getObjectVariable(rtrow,cand,'pt')
            etas[cand] = self.getObjectVariable(rtrow,cand,'eta')
            phis[cand] = self.getObjectVariable(rtrow,cand,'phi')
            p4s[cand] = self.getObjectVariable(rtrow,cand,'p4')
            charges[cand] = self.getObjectVariable(rtrow,cand,'charge')

        # require ++- or --+
        if abs(sum([charges[c] for c in leps]))!=1: return candidate

        # get the candidates
        hppCands = []
        for trio in itertools.permutations(leps,3):
            # require ++-/--+
            if charges[trio[0]]!=charges[trio[1]]: continue
            # require deltaR seperation of 0.02
            keep = True
            for i,j in itertools.combinations(range(3),2):
                dr = deltaR(etas[trio[i]],phis[trio[i]],etas[trio[j]],phis[trio[j]])
                m = (p4s[trio[i]]+p4s[trio[j]]).M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            if not keep: continue
            # require lead e/m pt > 20, if all taus, lead 2 pt>35
            ems = [cand for cand in trio if cand[0] in ['electrons','muons']]
            ts = [cand for cand in trio if cand[0] in ['taus']]
            if len(ems)>0:
                pts_ems = [pts[cand] for cand in ems]
                if max(pts_ems)<20.: continue
            else:
                pts_ts = [pts[cand] for cand in ts]
                if sorted(pts_ts)[-2]<35.: continue
            hppCands += [trio]
        if not hppCands: return candidate
        hppCand = hppCands[0][:2]
        hmCand = hppCands[0][2]

        candidate['hpp1'] = hppCand[0] if pts[hppCand[0]]>pts[hppCand[1]] else hppCand[1]
        candidate['hpp2'] = hppCand[1] if pts[hppCand[0]]>pts[hppCand[1]] else hppCand[0]
        candidate['hm1'] = hmCand

        # add jet
        jets = self.getCands(rtrow, 'jets', lambda rtrow,cand: self.getObjectVariable(rtrow,cand,'isLoose')>0.5)
        if len(jets)==1:
            candidate['leadJet'] = jets[0]
            candidate['subleadJet'] = ('jets',-1)
        if len(jets)>1:
            candidate['leadJet'] = jets[0]
            candidate['subleadJet'] = jets[1]
        else:
            candidate['leadJet'] = ('jets',-1)
            candidate['subleadJet'] = ('jets',-1)

        # do the z alternative combination
        massDiffs = {}
        for zpair in itertools.combinations(pts.keys(),2):
            if zpair[0][0]!=zpair[1][0]: continue # SF
            if charges[zpair[0]]==charges[zpair[1]]: continue # OS
            zp4 = p4s[zpair[0]] + p4s[zpair[1]]
            zmass = zp4.M()
            massDiffs[zpair] = abs(zmass-ZMASS)

        if not massDiffs:
            candidate['z1'] = ('',-1)
            candidate['z2'] = ('',-1)
            candidate['w1'] = ('',-1)
            return candidate

        # sort by closest z
        bestZ = sorted(massDiffs.items(), key=operator.itemgetter(1))[0][0]

        # now get the highest pt w
        zpts = {}
        zpts[bestZ[0]] = pts.pop(bestZ[0])
        zpts[bestZ[1]] = pts.pop(bestZ[1])
        bestW = sorted(pts.items(), key=operator.itemgetter(1))[-1][0]

        # and sort pt of Z
        z = sorted(zpts.items(), key=operator.itemgetter(1))
        z1 = z[1][0]
        z2 = z[0][0]

        candidate['z1'] = z1
        candidate['z2'] = z2
        candidate['w1'] = bestW

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

    def mediumFakeRate(self,rtrow,cand):
        return self.fakeRates.getFakeRate(rtrow,cand,'WZMedium','WZLoose')

    def tightFakeRate(self,rtrow,cand):
        return self.fakeRates.getFakeRate(rtrow,cand,'WZTight','WZLoose')

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
        for coll in ['electrons','muons','taus']:
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

    def numCentralJets(self,rtrow,cands,mode,pt):
        if cands['leadJet'][1]<0: return -1
        if cands['subleadJet'][1]<0: return -1
        eta1 = self.getObjectVariable(rtrow,cands['leadJet'],'eta')
        eta2 = self.getObjectVariable(rtrow,cands['subleadJet'],'eta')
        mineta = min(eta1,eta2)
        maxeta = max(eta1,eta2)
        return len(
            self.getCands(
                rtrow,
                'jets',
                lambda rtrow,cand: self.getObjectVariable(rtrow,cand,mode)>0.5
                                   and self.getObjectVariable(rtrow,cand,'pt')>pt
                                   and self.getObjectVariable(rtrow,cand,'eta')>mineta
                                   and self.getObjectVariable(rtrow,cand,'eta')<maxeta
            )
        )
    
    def zeppenfeld(self,rtrow,cands,*probeCands):
        if cands['leadJet'][1]<0: return -10.
        if cands['subleadJet'][1]<0: return -10.
        eta1 = self.getObjectVariable(rtrow,cands['leadJet'],'eta')
        eta2 = self.getObjectVariable(rtrow,cands['subleadJet'],'eta')
        meaneta = (eta1+eta2)/2
        if len(probeCands)>1:
            eta = self.getCompositeVariable(rtrow,'eta',*probeCands)
        else:
            eta = self.getObjectVariable(rtrow,probeCands[0],'eta')
        return eta-meaneta

    ######################
    ### channel string ###
    ######################
    def getChannelString(self,rtrow,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['hpp1','hpp2','hm1']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    def getWZChannelString(self,rtrow,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z1','z2','w1']:
            if cands[c][1]<0: return 'aaa'
            chanString += self.getCollectionString(cands[c])
        return chanString

    def getGenChannelString(self,rtrow,cands):
        '''Get the gen h++ channel'''
        chanString = ''
        pdgMap = {
            11: 'e',
            13: 'm',
            15: 't',
        }
        if 'HPlusPlusHMinusMinusHTo4L' in self.fileNames[0]: # h++h-- signal sample
            for s in [1,-1]:
                h = -1*s*9900041                     # h++ in pythia8
                for l1 in [s*11, s*13, s*15]:        # lepton 1
                    for l2 in [s*11, s*13, s*15]:    # lepton 2
                        if abs(l2)<abs(l1): continue # skip double counting
                        hasDecay = self.findDecay(rtrow,h,l1,l2)
                        if hasDecay:
                            chanString += pdgMap[abs(l1)]
                            chanString += pdgMap[abs(l2)]
        else:
            chanString = 'a'
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,rtrow,cands):
        return len(self.getPassingCands(rtrow,'Loose'))>=3

    def vetoFourth(self,rtrow,cands):
        return len(self.getPassingCands(rtrow,'Medium'))<=3

    def trigger(self,rtrow,cands):
        # accept MC, check trigger for data
        if rtrow.isData<0.5: return True
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
        candList = [cands[c] for c in ['hpp1','hpp2','hm1']]
        numTaus = [c[0] for c in candList].count('taus')
        if numTaus<2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12']
        elif numTaus==2:
            triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','DoublePFTau35']
        else:
            triggerList = ['DoublePFTau35']
        return self.triggerScales.getDataEfficiency(rtrow,triggerList,candList)







