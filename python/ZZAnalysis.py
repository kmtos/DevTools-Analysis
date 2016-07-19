#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR
from leptonId import passWZLoose, passWZMedium, passWZTight, passHppLoose, passHppMedium, passHppTight, passHZZLoose, passHZZTight

from Candidates import *

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("ZZAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


class ZZAnalysis(AnalysisBase):
    '''
    ZZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','zzlTree.root')
        outputTreeName = kwargs.pop('outputTreeName','ZZTree')
        self.preselection = 'muons_count+electrons_count>3'
        super(ZZAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.fourLoose,'fourLooseLeptons')
        self.cutTree.add(self.trigger,'trigger')

        # setup analysis tree

        # chan string
        self.tree.add(self.getChannelString, 'channel', ['C',5])

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
        else:
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVLPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.IsoMu22Pass(), 'pass_IsoMu22', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu22Pass(), 'pass_IsoTkMu22', 'I')
            self.tree.add(lambda cands: self.event.Mu45_eta2p1Pass(), 'pass_Mu45_eta2p1', 'I')
            self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.Ele25_eta2p1_WPTight_GsfPass(), 'pass_Ele25_eta2p1_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele27_eta2p1_WPLoose_GsfPass(), 'pass_Ele27_eta2p1_WPLoose_Gsf', 'I')
            self.tree.add(lambda cands: self.event.Ele45_WPLoose_GsfPass(), 'pass_Ele45_WPLoose_Gsf', 'I')
        self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
        self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
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
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z11'],cands['z12'],cands['z21'],cands['z22']), '4l_zeppenfeld','F')

        # z1 leptons
        self.addDiLepton('z1')
        self.addCompositeMet('z1met')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z11'],cands['z12']), 'z1_zeppenfeld','F')
        self.addLepton('z11')
        self.tree.add(lambda cands: self.passMedium(cands['z11']), 'z11_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z11']), 'z11_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z11']), 'z11_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z11']), 'z11_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z11']), 'z11_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z11']), 'z11_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z11']), 'z11_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z11']), 'z11_zeppenfeld','F')
        self.addLepton('z12')
        self.tree.add(lambda cands: self.passMedium(cands['z12']), 'z12_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z12']), 'z12_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z12']), 'z12_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z12']), 'z12_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z12']), 'z12_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z12']), 'z12_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z12']), 'z12_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z12']), 'z12_zeppenfeld','F')

        # z2 leptons
        self.addDiLepton('z2')
        self.addCompositeMet('z2met')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z21'],cands['z22']), 'z2_zeppenfeld','F')
        self.addLepton('z21')
        self.tree.add(lambda cands: self.passMedium(cands['z21']), 'z21_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z21']), 'z21_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z21']), 'z21_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z21']), 'z21_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z21']), 'z21_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z21']), 'z21_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z21']), 'z21_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z21']), 'z21_zeppenfeld','F')
        self.addLepton('z22')
        self.tree.add(lambda cands: self.passMedium(cands['z22']), 'z22_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z22']), 'z22_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z22']), 'z22_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z22']), 'z22_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z22']), 'z22_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z22']), 'z22_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z22']), 'z22_tightFakeRate', 'F')
        #self.tree.add(lambda cands: self.zeppenfeld(cands,cands['z22']), 'z22_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('z21_z11')
        self.addDiLepton('z21_z12')
        self.addDiLepton('z22_z11')
        self.addDiLepton('z22_z12')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z11' : None,
            'z12' : None,
            'z21' : None,
            'z22' : None,
            'z1' : None,
            'z1met' : None,
            'z2' : None,
            'z2met' : None,
            '4l' : None,
            '4lmet' : None,
            'z21_z11' : None,
            'z21_z12' : None,
            'z22_z11' : None,
            'z22_z12' : None,
            #'leadJet' : (),
            #'subleadJet' : (),
            'met': self.pfmet,
            'cleanJets' : [],
        }

        # get leptons
        leps = self.getPassingCands('Loose')
        medLeps = self.getPassingCands('Medium')
        if len(leps)<4: return candidate # need at least 4 leptons


        # get the candidates
        zzCands = []
        for quad in itertools.permutations(leps,4):
            # require +-+-
            if quad[0].charge()+quad[1].charge()!=0: continue
            if quad[2].charge()+quad[3].charge()!=0: continue
            # require same type
            if quad[0].collName!=quad[1].collName: continue
            if quad[2].collName!=quad[3].collName: continue
            # require deltaR seperation of 0.02 m(ll)>12
            keep = True
            for i,j in itertools.combinations(range(4),2):
                dicand = DiCandidate(quad[i],quad[j])
                dr = dicand.deltaR()
                m = dicand.M()
                if dr<0.02: keep = False
                if m<12.: keep = False
            if not keep: continue
            # require lead e/m pt > 25
            ems = [cand for cand in quad if cand.collName in ['electrons','muons']]
            if len(ems)>0:
                pts_ems = sorted([cand.pt() for cand in ems])
                if pts_ems[-1]<30.: continue
                if len(ems)>1:
                    if pts_ems[-2]<20.: continue
            # its a good candidate
            zzCands += [quad]
        if not zzCands: return candidate

        # sort by best z then st
        bestZ = 0
        highestSt = 0
        bestCand = []
        for quad in zzCands:
            z = DiCandidate(quad[0],quad[1])
            zmass = z.M()
            st = sum([cand.pt() for cand in quad[2:]])
            if abs(ZMASS-zmass)<(ZMASS-bestZ):
                bestZ = zmass
                highestSt = st
                bestCand = quad
            elif abs(ZMASS-zmass)==(ZMASS-bestZ):
                bestZ = zmass
                highestSt = st
                bestCand = quad

        z11 = bestCand[0] if bestCand[0].pt()>bestCand[1].pt() else bestCand[1]
        z12 = bestCand[1] if bestCand[0].pt()>bestCand[1].pt() else bestCand[0]
        z21 = bestCand[2] if bestCand[2].pt()>bestCand[3].pt() else bestCand[3]
        z22 = bestCand[3] if bestCand[2].pt()>bestCand[3].pt() else bestCand[2]

        candidate['z11'] = z11
        candidate['z12'] = z12
        candidate['z21'] = z21
        candidate['z22'] = z22
        candidate['z1'] = DiCandidate(z11,z12)
        candidate['z1met'] = MetCompositeCandidate(self.pfmet,z11,z12)
        candidate['z2'] = DiCandidate(z21,z22)
        candidate['z2met'] = MetCompositeCandidate(self.pfmet,z21,z22)
        candidate['4l'] = CompositeCandidate(z11,z12,z21,z22)
        candidate['4lmet'] = MetCompositeCandidate(self.pfmet,z11,z12,z21,z22)
        candidate['z21_z11'] = DiCandidate(z21,z11)
        candidate['z21_z12'] = DiCandidate(z21,z12)
        candidate['z22_z11'] = DiCandidate(z22,z11)
        candidate['z22_z12'] = DiCandidate(z22,z12)

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps,0.4)

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
            #return self.leptonScales.getScale('None',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedVeto',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDVeto',cand)
            #return self.leptonScales.getScale('None',cand)
        else:
            return 1.

    def mediumScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
            #return self.leptonScales.getScale('HZZTight',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedMedium',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDMedium',cand)
            #return self.leptonScales.getScale('HZZTight',cand)
        else:
            return 1.

    def tightScale(self,cand):
        if cand.collName=='muons':
            return self.leptonScales.getScale('MediumIDTightIso',cand)
            #return self.leptonScales.getScale('HZZTight',cand)
        elif cand.collName=='electrons':
            return self.leptonScales.getScale('CutbasedTight',cand) if self.version=='76X' else self.leptonScales.getScale('CutBasedIDTight',cand)
            #return self.leptonScales.getScale('HZZTight',cand)
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
        for coll in [self.electrons,self.muons]:
            cands += self.getCands(coll,passMode)
        return cands

    def numJets(self,cleanJets,mode,pt):
        jetColl = self.getCands(
            cleanJets,
            lambda cand: getattr(cand,mode)()>0.5 and cand.pt()>pt
        )
        return len(jetColl)


    ######################
    ### channel string ###
    ######################
    def getChannelString(self,cands):
        '''Get the channel string'''
        chanString = ''
        for c in ['z11','z12','z21','z22']:
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
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    'IsoMu22',
                    'IsoTkMu22',
                    'Mu45_eta2p1',
                    'Mu50',
                ],
                'SingleElectron' : [
                    'Ele25_eta2p1_WPTight_Gsf',
                    'Ele27_WPTight_Gsf',
                    'Ele27_eta2p1_WPLoose_Gsf',
                    'Ele45_WPLoose_Gsf',
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
        candList = [cands[c] for c in ['z11','z12','z21','z22']]
        triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['SingleMuSoup','SingleEleSoup','Mu17Mu8','Ele23Ele12']
        return self.triggerScales.getDataEfficiency(triggerList,candList)








def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('zz'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='zzTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    zzAnalysis = ZZAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='ZZTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       zzAnalysis.analyze()
       zzAnalysis.finish()
    except KeyboardInterrupt:
       zzAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
