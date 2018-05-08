#!/usr/bin/env python
import argparse
import logging
import sys

from DevTools.Analyzer.utilities import getTestFiles
from AnalysisBase import AnalysisBase
from utilities import ZMASS, deltaPhi, deltaR

from Candidates import *
import KinematicFitter

import sys
import itertools
import operator

import ROOT

logger = logging.getLogger("MuMuTauTauAnalysis")
logging.basicConfig(level=logging.INFO, stream=sys.stderr,format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def load_events(h,a):
    events = []
    #with open('events_{h}_{a}_ktos.txt'.format(h=h,a=a)) as f:
    with open('h{h}a{a}_Events_in_Kyle_Not_Devin.txt'.format(h=h,a=a)) as f:
        for l in f.readlines():
            events += [l.strip()]
    return events

class MuMuTauTauAnalysis(AnalysisBase):
    '''
    MuMuTauTau analysis
    '''

    def __init__(self,**kwargs):
        outputFileName = kwargs.pop('outputFileName','muMuTauTauTree.root')
        outputTreeName = kwargs.pop('outputTreeName','MuMuTauTauTree')
        super(MuMuTauTauAnalysis, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,**kwargs)

        self._kinfit = None

        # setup cut tree
        self.cutTree.add(self.metFilter,'metFilter')
        self.cutTree.add(self.trigger,'trigger')

        # open a file of events
        self.events = []
        if 'SUSYGluGluToHToAA_AToMuMu_AToTauTau' in self.fileNames[0]:
            for h in [125,300,750]:
                if 'M-{h}_'.format(h=h) not in self.fileNames[0] and h!=125: continue
                for a in ['3p6',4,5,6,7,9,11,13,15,17,19,21]:
                    if 'M-{a}_'.format(a=a) not in self.fileNames[0]: continue
                    try:
                        self.events = load_events(h,a)
                    except:
                        logging.warning('failed to load events {h} {a}'.format(h=h,a=a))

        # setup analysis tree

        # event counts
        self.tree.add(self.getGenChannel, 'genChannel', ['C',5])
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isLoose',20), 'numJetsLoose20', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'isTight',20), 'numJetsTight20', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2L',20), 'numBjetsLoose20', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2M',20), 'numBjetsMedium20', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJets'],'passCSVv2T',20), 'numBjetsTight20', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJetsDR08'],'isLoose',20), 'numJetsLoose20DR08', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJetsDR08'],'isTight',20), 'numJetsTight20DR08', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJetsDR08'],'passCSVv2L',20), 'numBjetsLoose20DR08', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJetsDR08'],'passCSVv2M',20), 'numBjetsMedium20DR08', 'I')
        self.tree.add(lambda cands: self.numJets(cands['cleanJetsDR08'],'passCSVv2T',20), 'numBjetsTight20DR08', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passLoose)), 'numLooseElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.electrons,self.passTight)), 'numTightElectrons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passLoose)), 'numLooseMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.muons,self.passTight)), 'numTightMuons', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passLoose)), 'numLooseTaus', 'I')
        self.tree.add(lambda cands: len(self.getCands(self.taus,self.passTight)), 'numTightTaus', 'I')

        # trigger
        if self.version=='76X':
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
        else:
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
            self.tree.add(lambda cands: self.event.Mu17_Mu8_SameSign_DZPass(), 'pass_Mu17_Mu8_SameSign_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu20_Mu10_SameSign_DZPass(), 'pass_Mu20_Mu10_SameSign_DZ', 'I')

        self.tree.add(lambda cands: self.triggerEfficiency(cands)[0], 'triggerEfficiency', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[1], 'triggerEfficiencyUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiency(cands)[2], 'triggerEfficiencyDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[0], 'triggerEfficiencyMC', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[1], 'triggerEfficiencyMCUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyMC(cands)[2], 'triggerEfficiencyMCDown', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[0], 'triggerEfficiencyData', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[1], 'triggerEfficiencyDataUp', 'F')
        self.tree.add(lambda cands: self.triggerEfficiencyData(cands)[2], 'triggerEfficiencyDataDown', 'F')

        # h
        self.addGenParticle('gh')
        self.addComposite('h')
        self.tree.add(lambda cands: self.kinfit(cands).getComposite('h').M(), 'h_massKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getFitStatus(), 'kinFitStatus', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).atLowerBound, 'kinFitAtLowerBound', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).atUpperBound, 'kinFitAtUpperBound', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).getFitX(), 'kinFitX', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getFitChi2(), 'kinFitChi2', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getRecoil().Px(), 'kinFitRecoilPx', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getRecoil().Py(), 'kinFitRecoilPy', 'F')
        self.addCompositeMet('hmet')
        self.tree.add(lambda cands: cands['hmet'].Mcat(2,3), 'hmet_mcat', 'F')

        # amm leptons
        self.addCompositeGenParticle('gamm','gam1','gam2')
        self.addDiLepton('amm')
        self.tree.add(lambda cands: self.kinfit(cands).getComposite('amm').M(), 'amm_massKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getCompositeStatus('amm'), 'amm_kinFitStatus', 'I')
        self.addCompositeMet('ammmet')
        self.addLepton('am1')
        self.addGenDeltaR('am1','gam1')
        self.addGenParticle('gam1')
        self.addDetailedMuon('am1')
        self.tree.add(lambda cands: self.kinfit(cands).getParticleStatus('m1'),    'am1_kinFitStatus', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m1').Pt(),     'am1_ptKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m1').Eta(),    'am1_etaKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m1').Phi(),    'am1_phiKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m1').Energy(), 'am1_energyKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m1').M(),      'am1_massKinFit', 'F')
        self.addLeptonMet('am1met')
        self.addLepton('am2')
        self.addGenDeltaR('am2','gam2')
        self.addGenParticle('gam2')
        self.addDetailedMuon('am2')
        self.tree.add(lambda cands: self.kinfit(cands).getParticleStatus('m2'),    'am2_kinFitStatus', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m2').Pt(),     'am2_ptKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m2').Eta(),    'am2_etaKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m2').Phi(),    'am2_phiKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m2').Energy(), 'am2_energyKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('m2').M(),      'am2_massKinFit', 'F')
        self.addLeptonMet('am2met')

        # att leptons
        self.addCompositeGenParticle('gatt','gat1','gat2')
        self.addDiLepton('att')
        self.tree.add(lambda cands: self.kinfit(cands).getComposite('att').M(), 'att_massKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getCompositeStatus('att'), 'att_kinFitStatus', 'I')
        self.addCompositeMet('attmet')
        self.tree.add(lambda cands: cands['attmet'].Mcat(0,1), 'attmet_mcat', 'F')
        self.addGenParticle('gat1',isTau=True)
        self.addLepton('atm')
        self.addGenDeltaR('atm','gatm')
        self.addDetailedMuon('atm')
        self.tree.add(lambda cands: self.kinfit(cands).getParticleStatus('tm'),    'atm_kinFitStatus', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('tm').Pt(),     'atm_ptKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('tm').Eta(),    'atm_etaKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('tm').Phi(),    'atm_phiKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('tm').Energy(), 'atm_energyKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('tm').M(),      'atm_massKinFit', 'F')
        self.addLeptonMet('atmmet')
        self.addGenParticle('gat2',isTau=True)
        self.addLepton('ath')
        self.addGenDeltaR('ath','gath')
        self.addDetailedTau('ath')
        self.tree.add(lambda cands: self.kinfit(cands).getParticleStatus('th'), 'ath_kinFitStatus', 'I')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('th').Pt(),     'ath_ptKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('th').Eta(),    'ath_etaKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('th').Phi(),    'ath_phiKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('th').Energy(), 'ath_energyKinFit', 'F')
        self.tree.add(lambda cands: self.kinfit(cands).getParticle('th').M(),      'ath_massKinFit', 'F')
        self.addLeptonMet('athmet')
        self.addJet('athjet')

        # met
        self.addMet('met')

        # other event
        self.tree.add(lambda cands: sum([x.pt() for x in cands['cleanJets']]), 'ht', 'F')

    ############################
    ### Additional functions ###
    ############################
    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        self.addCandVar(label,'matches_IsoMu24','matches_IsoMu24','I')
        self.addCandVar(label,'matches_IsoTkMu24','matches_IsoTkMu24','I')
        super(MuMuTauTauAnalysis,self).addDetailedMuon(label)

    def addGenParticle(self,label,isTau=False):
        self.tree.add(lambda cands: cands[label].pt()     if label in cands else 0, '{0}_pt'.format(label),     'F')
        self.tree.add(lambda cands: cands[label].eta()    if label in cands else 0, '{0}_eta'.format(label),    'F')
        self.tree.add(lambda cands: cands[label].phi()    if label in cands else 0, '{0}_phi'.format(label),    'F')
        self.tree.add(lambda cands: cands[label].energy() if label in cands else 0, '{0}_energy'.format(label), 'F')
        self.tree.add(lambda cands: cands[label].mass()   if label in cands else 0, '{0}_mass'.format(label),   'F')
        self.tree.add(lambda cands: cands[label].charge() if label in cands else 0, '{0}_charge'.format(label), 'F')
        self.tree.add(lambda cands: cands[label].pdgId()  if label in cands else 0, '{0}_pdgId'.format(label),  'I')
        if isTau:
            self.tree.add(lambda cands: self.getGenDecayMode(cands[label])  if label in cands else 0, '{0}_decayMode'.format(label),  'I')

    def addCompositeGenParticle(self,label,d1,d2):
        self.addGenParticle(label)
        self.tree.add(lambda cands: deltaR(cands[d1].eta(),cands[d1].phi(),cands[d2].eta(),cands[d2].phi()) if label in cands else 0, '{0}_deltaR'.format(label), 'F')

    def addGenDeltaR(self,label,gen):
        self.tree.add(lambda cands: deltaR(cands[label].eta(),cands[label].phi(),cands[gen].eta(),cands[gen].phi()) if gen in cands else 9, '{0}_genTruthDeltaR'.format(label), 'F')

    def passMuon(self,cand):
        if cand.pt()<3: return False
        if not cand.isPFMuon(): return False
        if not (cand.isGlobalMuon() or cand.isTrackerMuon()): return False
        if abs(cand.dxy())>=0.2: return False
        if abs(cand.dz())>=0.5: return False
        return True

    def passTau(self,cand):
        if cand.pt()<10: return False
        if abs(cand.dxy())>=0.2: return False
        if abs(cand.dz())>=0.5: return False
        if not cand.decayModeFinding(): return False
        return True

    def kinfit(self,cands):
        # return fitted constraints
        if self._kinfit: return self._kinfit

        # create constraints
        m1 = cands['am1']
        m2 = cands['am2']
        tm = cands['atm']
        th = cands['ath']
        m1p4 = KinematicFitter.Muon(       'm1', m1.pt(), m1.eta(), m1.phi(), m1.energy())
        m2p4 = KinematicFitter.Muon(       'm2', m2.pt(), m2.eta(), m2.phi(), m2.energy())
        tmp4 = KinematicFitter.MuonTau(    'tm', tm.pt(), tm.eta(), tm.phi(), tm.energy())
        thp4 = KinematicFitter.HadronicTau('th', th.pt(), th.eta(), th.phi(), th.energy(), th.decayMode())
        # TODO: add cov of muons/taus with appropriate error on energy
        m1p4.setErrors(0.01*m1p4.E(),0,0)
        m2p4.setErrors(0.01*m2p4.E(),0,0)

        met = cands['met']
        metp4 = KinematicFitter.MET('met', met.et(), met.phi(), met.cov00(), met.cov01(), met.cov01(), met.cov11())

        kinfit = KinematicFitter.KinematicFitter()
        kinfit.addParticle('m1',m1p4)
        kinfit.addParticle('m2',m2p4)
        kinfit.addParticle('tm',tmp4)
        kinfit.addParticle('th',thp4)
        kinfit.addParticle('met',metp4)

        kinfit.addComposite('amm','m1','m2')
        kinfit.addComposite('att','tm','th')
        kinfit.addComposite('h','m1','m2','tm','th')

        kinfit.addMassConstraint('amm','att',0)

        # perform fit with constraints
        kinfit.setMinimizationParticles('th')
        kinfit.fit()

        self._kinfit = kinfit

        return self._kinfit

    ############################
    ### select 4l candidates ###
    ############################
    def selectCandidates(self):
        # reset kinfit
        self._kinfit = None

        candidate = {
            'am1' : None,
            'am1met' : None,
            'am2' : None,
            'am2met' : None,
            'atm' : None,
            'atmmet' : None,
            'ath' : None,
            'athmet' : None,
            'amm' : None,
            'ammmet' : None,
            'att' : None,
            'attmet' : None,
            'h' : None,
            'hmet' : None,
            'met': self.pfmet,
            'athjet': Candidate(None),
            'cleanJets' : [],
            'cleanJetsNoTau' : [],
        }

        # get leptons
        muons = [m for m in self.muons if self.passMuon(m)]
        taus = [t for t in self.taus if self.passTau(t)]
        leps = muons+taus
        if len(muons)<3:
            self.report_failure('fails 3 muon requirement')
            return candidate
        if len(taus)<1: 
            self.report_failure('fails 1 tau requirement')
            return candidate


        # get the candidates
        hCand = []
        mmDeltaR = 999
        ttDeltaR = 999
        furthest = 0
        nperms = 0
        nvalid = 0
        for quad in itertools.permutations(leps,4):
            nperms += 1
            # require mmmt
            if not quad[0].__class__.__name__=='Muon': continue
            if not quad[1].__class__.__name__=='Muon': continue
            if not quad[2].__class__.__name__=='Muon': continue
            if not quad[3].__class__.__name__=='Tau': continue
            # trigger match
            matchTrigger = quad[0].matches_IsoMu24() or quad[0].matches_IsoTkMu24()
            if not matchTrigger: continue
            furthest = max([1,furthest])
            # charge OS
            if quad[0].charge()==quad[1].charge(): continue
            if quad[2].charge()==quad[3].charge(): continue
            furthest = max([2,furthest])
            nvalid += 1
            # require lead m pt>25
            if quad[0].pt()<25: continue
            furthest = max([3,furthest])
            # make composites
            amm = DiCandidate(quad[0],quad[1])
            att = DiCandidate(quad[2],quad[3])
            # reject events with mu1/mu2 near tm/th
            m1tm = DiCandidate(quad[0],quad[2])
            m1th = DiCandidate(quad[0],quad[3])
            m2tm = DiCandidate(quad[1],quad[2])
            m2th = DiCandidate(quad[1],quad[3])
            if m1tm.deltaR()<0.4: continue
            if m1th.deltaR()<0.8: continue
            if m2tm.deltaR()<0.4: continue
            if m2th.deltaR()<0.8: continue
            furthest = max([4,furthest])
            # choose best
            if not hCand: hCand = quad
            better = True
            if amm.deltaR()>mmDeltaR:
                better = False
            elif amm.deltaR()==mmDeltaR and att.deltaR()>ttDeltaR:
                better = False
            if better:
                hCand = quad
                mmDeltaR = amm.deltaR()
                ttDeltaR = att.deltaR()

        #print 'number perms: {}, number valid: {}'.format(nperms,nvalid)

        furthestMap = {
            0: 'topology',
            1: 'trigger',
            2: 'OS',
            3: 'lead pt',
            4: 'deltaR',
        }
                
        if not hCand:
            self.report_failure('no higgs candidate, furthest {}'.format(furthestMap[furthest]))
            return candidate

        am1 = hCand[0] if hCand[0].pt()>hCand[1].pt() else hCand[1]
        am2 = hCand[1] if hCand[0].pt()>hCand[1].pt() else hCand[0]
        atm = hCand[2]
        ath = hCand[3]

        amm = DiCandidate(am1,am2)
        if amm.M()>30:
            self.report_failure('selected higgs amm M>30')
            return candidate
        #if amm.deltaR()>1.5: return candidate

        candidate['am1'] = am1
        candidate['am1met'] = MetCompositeCandidate(self.pfmet,am1)
        candidate['am2'] = am2
        candidate['am2met'] = MetCompositeCandidate(self.pfmet,am2)
        candidate['atm'] = atm
        candidate['atmmet'] = MetCompositeCandidate(self.pfmet,atm)
        candidate['ath'] = ath
        candidate['athmet'] = MetCompositeCandidate(self.pfmet,ath)
        candidate['amm'] = DiCandidate(am1,am2)
        candidate['ammmet'] = MetCompositeCandidate(self.pfmet,am1,am2)
        candidate['att'] = DiCandidate(atm,ath)
        candidate['attmet'] = MetCompositeCandidate(self.pfmet,atm,ath)
        candidate['h'] = CompositeCandidate(am1,am2,atm,ath)
        candidate['hmet'] = MetCompositeCandidate(self.pfmet,am1,am2,atm,ath)

        # clean the jets
        candidate['cleanJets'] = self.cleanCands(self.jets,[am1,am2,atm,ath],0.4)
        candidate['cleanJetsDR08'] = self.cleanCands(candidate['cleanJets'],[ath],0.8)

        # match jet to tau
        dr = 999
        j = None
        for jet in self.jets:
            jt = DiCandidate(jet,ath)
            if jt.deltaR()<dr:
                j = jet
                dr = jt.deltaR()
        if j:
            candidate['athjet'] = j

        candidate.update(self.getGenCandidates())

        return candidate

    def getGenCandidates(self):
        if 'SUSYGluGluToHToAA_AToMuMu_AToTauTau' not in self.fileNames[0]:
            return {}
        gmuons = [g for g in self.gen if abs(g.pdgId())==13]
        gtaus = [g for g in self.gen if abs(g.pdgId())==15]
        gas = [g for g in self.gen if abs(g.pdgId())==36]
        ghs = [g for g in self.gen if abs(g.pdgId()) in [25,35]]

        gmfromas = [g for g in gmuons if g.mother_1()==36 or g.mother_2()==36]
        gtfromas = [g for g in gtaus if g.mother_1()==36 or g.mother_2()==36]
        lastgh = [g for g in ghs if g.daughter_1()==36 and g.daughter_2()==36]

        if len(lastgh)<1:
            logging.warning('No h cand found')

        h = lastgh[0]
        amm = gas[0] if abs(gas[0].daughter_1())==13 else gas[1]
        att = gas[1] if abs(gas[0].daughter_1())==13 else gas[0]
        am1 = gmfromas[0] if gmfromas[0].pt()>gmfromas[1].pt() else gmfromas[1]
        am2 = gmfromas[1] if gmfromas[0].pt()>gmfromas[1].pt() else gmfromas[0]
        at1 = gtfromas[0] if gtfromas[0].pt()>gtfromas[1].pt() else gtfromas[1]
        at2 = gtfromas[1] if gtfromas[0].pt()>gtfromas[1].pt() else gtfromas[0]

        cands = {
            'gh': h,
            'gamm': amm,
            'gam1': am1,
            'gam2': am2,
            'gatt': att,
            'gat1': at1,
            'gat2': at2,
        }

        if self.getGenDecayMode(at1)==-13 and self.getGenDecayMode(at2)==99:
            cands['gatm'] = at1
            cands['gath'] = at2
        elif self.getGenDecayMode(at2)==-13 and self.getGenDecayMode(at1)==99:
            cands['gatm'] = at2
            cands['gath'] = at1
        else:
            if self.getGenDecayMode(at2)==-13:
                cands['gatm'] = at2
            if self.getGenDecayMode(at2)==99:
                cands['gath'] = at2
            if self.getGenDecayMode(at1)==-13:
                cands['gatm'] = at1
            if self.getGenDecayMode(at1)==99:
                cands['gath'] = at1


        return cands

    def getGenChannel(self,cands):
        if 'gh' not in cands: return 'xxxx'
        dm1 = self.getGenDecayMode(cands['gat1'])
        dm2 = self.getGenDecayMode(cands['gat2'])
        channel = 'mm'
        dmMap = {
            -11: 'e',
            -13: 'm',
            -99: 'x',
            99 : 'h',
            0  : 'x',
        }
        channel += dmMap[dm1]
        channel += dmMap[dm2]

        return channel

    # find tau decay modes
    def getGenDecayMode(self,t):
        gtaus = [g for g in self.gen if abs(g.pdgId())==15]
        # check for decays to taus
        if abs(t.daughter_1())==15 or abs(t.daughter_2())==15:
            # find the last daughter that doesnt decay to a tau
            oktaus = [g for g in gtaus if g.pdgId()==t.pdgId()]
            oktaus = [g for g in oktaus if abs(g.daughter_1())!=15 and abs(g.daughter_2())!=15]
            if len(oktaus)==0: return -99
            return self.getGenDecayMode(oktaus[-1])
        # leptons
        elif abs(t.daughter_1())==13 or abs(t.daughter_2())==13:
            return -13
        elif abs(t.daughter_1())==11 or abs(t.daughter_2())==11:
            return -11
        # hadronic
        # not enough info to determine which dm, only have num particles here
        else:
            return 99


    ###########################
    ### analysis selections ###
    ###########################
    def trigger(self,cands):
        isData = self.event.isData()>0.5
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
                    'IsoMu24',
                    'IsoTkMu24',
                    #'Mu45_eta2p1',
                    #'Mu50',
                ],
            }
        # the order here defines the heirarchy
        # first dataset, any trigger passes
        # second dataset, if a trigger in the first dataset is found, reject event
        # so forth
        datasets = [
            'SingleMuon',
        ]
        result = self.checkTrigger(*datasets,**triggerNames)
        if not result: self.report_failure('fails trigger')
        return result

    def triggerEfficiencyMC(self,cands):
        return self.triggerEfficiency(cands,mode='mc')

    def triggerEfficiencyData(self,cands):
        return self.triggerEfficiency(cands,mode='data')

    def triggerEfficiency(self,cands,mode='ratio'):
        candList = [cands[c] for c in ['am1','am2','atm','ath']]
        triggerList = ['IsoMu20_OR_IsoTkMu20'] if self.version=='76X' else ['IsoMu24_OR_IsoTkMu24']
        if mode=='data':
            return self.triggerScales.getDataEfficiency(triggerList,candList,doError=True)
        elif mode=='mc':
            return self.triggerScales.getMCEfficiency(triggerList,candList,doError=True)
        elif mode=='ratio':
            return self.triggerScales.getRatio(triggerList,candList,doError=True)







def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('haa',version='80XMuMuTauTau'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='muMuTauTauTree.root', help='Output file')
    parser.add_argument('--shift', type=str, default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    muMuTauTauAnalysis = MuMuTauTauAnalysis(
        outputFileName=args.outputFile,
        outputTreeName='MuMuTauTauTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='MiniTree',
        inputLumiName='LumiTree',
        inputTreeDirectory='miniTree',
        shift = args.shift,
    )

    try:
       muMuTauTauAnalysis.analyze()
       muMuTauTauAnalysis.finish()
    except KeyboardInterrupt:
       muMuTauTauAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
