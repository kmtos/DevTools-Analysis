#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *

import itertools
import operator

import ROOT

logger = logging.getLogger("WZAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class WZAnalysis(AnalysisBase):
    '''
    WZ analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','wzTree.root')
        outputTreeName = kwargs.pop('outputTreeName','WZTree')
        super(WZAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        # setup cut tree
        self.cutTree.add(self.threeLoose,'threeLooseLeptons')
        self.cutTree.add(self.vetoFourth,'noFourthMediumLepton')
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

        # trigger
        self.addTriggers()
        self.tree.add(self.triggerEfficiency, 'triggerEfficiency', 'F')
        self.tree.add(self.triggerEfficiencyMC, 'triggerEfficiencyMC', 'F')
        self.tree.add(self.triggerEfficiencyData, 'triggerEfficiencyData', 'F')

        # vbf
        self.addJet('leadJet')
        self.addJet('subleadJet')
        self.addDiJet('dijet')
        self.tree.add(lambda cands: self.numCentralJets(cands['leadJet'],cands['subleadJet'],cands['cleanJets'],'isLoose',30), 'dijet_numCentralJetsLoose30', 'I')
        self.tree.add(lambda cands: self.numCentralJets(cands['leadJet'],cands['subleadJet'],cands['cleanJets'],'isTight',30), 'dijet_numCentralJetsTight30', 'I')

        # 3 lepton
        self.addComposite('3l')
        self.addCompositeMet('3lmet')
        self.tree.add(lambda cands: self.zeppenfeld(cands['leadJet'],cands['subleadJet'],cands['z1'],cands['z2'],cands['w1']), '3l_zeppenfeld','F')

        # z leptons
        self.addDiLepton('z')
        self.tree.add(lambda cands: self.zeppenfeld(cands['leadJet'],cands['subleadJet'],cands['z1'],cands['z2']), 'z_zeppenfeld','F')
        self.addLepton('z1')
        self.tree.add(lambda cands: self.passMedium(cands['z1']), 'z1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z1']), 'z1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z1'])[0], 'z1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z1'])[0], 'z1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z1'])[0], 'z1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z1'])[0], 'z1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z1'])[0], 'z1_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.zeppenfeld(cands['leadJet'],cands['subleadJet'],cands['z1']), 'z1_zeppenfeld','F')
        self.addLepton('z2')
        self.tree.add(lambda cands: self.passMedium(cands['z2']), 'z2_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['z2']), 'z2_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['z2'])[0], 'z2_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['z2'])[0], 'z2_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['z2'])[0], 'z2_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['z2'])[0], 'z2_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['z2'])[0], 'z2_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.zeppenfeld(cands['leadJet'],cands['subleadJet'],cands['z2']), 'z2_zeppenfeld','F')

        # w lepton
        self.addLeptonMet('w')
        self.addLepton('w1')
        self.tree.add(lambda cands: self.passMedium(cands['w1']), 'w1_passMedium', 'I')
        self.tree.add(lambda cands: self.passTight(cands['w1']), 'w1_passTight', 'I')
        self.tree.add(lambda cands: self.looseScale(cands['w1'])[0], 'w1_looseScale', 'F')
        self.tree.add(lambda cands: self.mediumScale(cands['w1'])[0], 'w1_mediumScale', 'F')
        self.tree.add(lambda cands: self.tightScale(cands['w1'])[0], 'w1_tightScale', 'F')
        self.tree.add(lambda cands: self.mediumFakeRate(cands['w1'])[0], 'w1_mediumFakeRate', 'F')
        self.tree.add(lambda cands: self.tightFakeRate(cands['w1'])[0], 'w1_tightFakeRate', 'F')
        self.tree.add(lambda cands: self.zeppenfeld(cands['leadJet'],cands['subleadJet'],cands['w1']), 'w1_zeppenfeld','F')

        # wrong combination
        self.addDiLepton('w1_z1')
        self.addDiLepton('w1_z2')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### select WZ candidates ###
    ############################
    def selectCandidates(self):
        candidate = {
            'z1' : None,
            'z2' : None,
            'w1' : None,
            'z' : None,
            'w' : None,
            '3l': None,
            '3lmet': None,
            'w1_z1': None,
            'w1_z2': None,
            'met': self.pfmet,
            'cleanJets': [],
            'leadJet': Candidate(None),
            'subleadJet': Candidate(None),
            'dijet': DiCandidate(Candidate(None),Candidate(None)),
        }

        #TODO: Fix low efficiency

        # get leptons
        leps = self.getPassingCands('Loose',self.electrons,self.muons)
        if len(leps)<3: return candidate # need at least 3 leptons

        # get invariant masses
        bestWZ = ()
        bestMassdiff = 99999
        bestWpt = 0
        for wz in itertools.permutations(leps,3):
            # Z selection
            if wz[0].collName!=wz[1].collName: continue # SF
            if wz[0].charge()==wz[1].charge(): continue # OS
            z = DiCandidate(*wz[:2])
            if z.deltaR()<0.02: continue
            # w to deltaR
            zw1 = DiCandidate(wz[0],wz[2])
            zw2 = DiCandidate(wz[1],wz[2])
            if zw1.deltaR()<0.02: continue
            if zw2.deltaR()<0.02: continue
            # choose best
            massdiff = abs(z.M()-ZMASS)
            wpt = wz[2].pt()
            if massdiff<bestMassdiff:
                bestWZ = wz
                bestMassdiff = massdiff
                bestWpt = wpt
            elif massdiff==bestMassdiff and wpt>bestWpt:
                bestWZ = wz
                bestMassdiff = massdiff
                bestWpt = wpt


        if not bestWZ: return candidate # need a wz candidate

        z1 = bestWZ[0] if bestWZ[0].pt()>bestWZ[1].pt() else bestWZ[1]
        z2 = bestWZ[1] if bestWZ[0].pt()>bestWZ[1].pt() else bestWZ[0]
        w1 = bestWZ[2]

        medLeps = self.getPassingCands('Medium',self.electrons,self.muons)
        remainingLeps = self.cleanCands(medLeps,[z1,z2,w1],0.02)

        if remainingLeps:
            #print '{0}:{1}:{2}'.format(self.event.run(),self.event.lumi(),self.event.event())
            #print 'Veto based on remaining Leps', len(remainingLeps)
            return candidate

        candidate['z1'] = z1
        candidate['z2'] = z2
        candidate['w1'] = w1
        candidate['z'] = DiCandidate(z1,z2)
        candidate['w'] = MetCompositeCandidate(self.pfmet,w1)
        candidate['3l'] = CompositeCandidate(z1,z2,w1)
        candidate['3lmet'] = MetCompositeCandidate(self.pfmet,z1,z2,w1)
        candidate['w1_z1'] = DiCandidate(w1,z1)
        candidate['w1_z2'] = DiCandidate(w1,z2)

        candidate['cleanJets'] = self.cleanCands(self.jets,medLeps+[z1,z2,w1],0.4)
        if len(candidate['cleanJets'])>0: candidate['leadJet'] = candidate['cleanJets'][0]
        if len(candidate['cleanJets'])>1: candidate['subleadJet'] = candidate['cleanJets'][1]
        if len(candidate['cleanJets'])>1: candidate['dijet'] = DiCandidate(candidate['leadJet'],candidate['subleadJet'])

        return candidate

    def numCentralJets(self,leadJet,subleadJet,cleanJets,mode,pt):
        if not leadJet: return -1
        if not subleadJet: return -1
        eta1 = leadJet.eta()
        eta2 = subleadJet.eta()
        mineta = min(eta1,eta2)
        maxeta = max(eta1,eta2)
        return len(
            self.getCands(
                cleanJets,
                lambda cand: getattr(cand,mode)()>0.5
                             and cand.pt()>pt
                             and cand.eta()>mineta
                             and cand.eta()<maxeta
            )
        )
    
    def zeppenfeld(self,leadJet,subleadJet,*probeCands):
        if not leadJet: return -10.
        if not subleadJet: return -10.
        eta1 = leadJet.eta()
        eta2 = subleadJet.eta()
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
        for c in ['z1','z2','w1']:
            chanString += self.getCollectionString(cands[c])
        return chanString

    ###########################
    ### analysis selections ###
    ###########################
    def threeLoose(self,cands):
        return len(self.getPassingCands('Loose',self.electrons,self.muons))>=3

    def vetoFourth(self,cands):
        return len(self.getPassingCands('Medium',self.electrons,self.muons))<=3

    def trigger(self,cands):
        # accept MC, check trigger for data
        #if self.event.isData()<0.5: return True
        # use trigger for MC and data
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
            }
        else:
            triggerNames = {
                'DoubleMuon'     : [
                    #'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL',
                    #'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL',
                    'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
                    'Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
                    'DiMu9_Ele9_CaloIdL_TrackIdL',
                    'TripleMu_12_10_5',
                ],
                'DoubleEG'       : [
                    'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                    'Ele16_Ele12_Ele8_CaloIdL_TrackIdL',
                    'Mu8_DiEle12_CaloIdL_TrackIdL',
                ],
                'MuonEG'         : [
                    'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
                    'Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL',
                ],
                'SingleMuon'     : [
                    #'IsoMu22',
                    #'IsoTkMu22',
                    'IsoMu24',
                    'IsoTkMu24',
                    #'Mu45_eta2p1',
                    #'Mu50',
                ],
                'SingleElectron' : [
                    'Ele25_eta2p1_WPTight_Gsf',
                    'Ele27_WPTight_Gsf',
                    #'Ele27_eta2p1_WPLoose_Gsf',
                    #'Ele45_WPLoose_Gsf',
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
        ]
        return self.checkTrigger(*datasets,**triggerNames)

    def triggerEfficiency(self,cands):
        return self.triggerScales.getDataEfficiency(triggerList,candList)

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['z1','z2','w1']]
        triggerList = ['IsoMu20_OR_IsoTkMu20','Ele23_WPLoose','Mu17_Mu8','Ele17_Ele12'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24','Ele27Tight','Mu17Mu8','Ele23Ele12']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList)








def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('wz'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='wzTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    wzAnalysis = WZAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='WZTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       wzAnalysis.analyze()
       wzAnalysis.finish()
    except KeyboardInterrupt:
       wzAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
