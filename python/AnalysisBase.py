# AnalysisBase.py

import logging
import os
import sys
import math
import time

sys.argv.append('-b')
import ROOT
sys.argv.pop()
from array import array

from CutTree import CutTree
from AnalysisTree import AnalysisTree

from PileupWeights import PileupWeights
from FakeRates import FakeRates
from LeptonScales import LeptonScales
from TriggerScales import TriggerScales
from TriggerPrescales import TriggerPrescales
from ZptGenWeight import ZptGenWeight

from utilities import deltaR, deltaPhi
from DevTools.Utilities.utilities import getCMSSWVersion
from Candidates import *

try:
    from progressbar import ProgressBar, ETA, Percentage, Bar, SimpleProgress
    hasProgress = False
except:
    hasProgress = False

class AnalysisBase(object):
    '''
    Analysis Tree
    '''

    def __init__(self,**kwargs):
        inputFileNames = kwargs.pop('inputFileNames',[])
        inputTreeDirectory = kwargs.pop('inputTreeDirectory','miniTree')
        inputTreeName = kwargs.pop('inputTreeName','MiniTree')
        inputLumiName = kwargs.pop('inputTreeName','LumiTree')
        outputFileName = kwargs.pop('outputFileName','analysisTree.root')
        outputTreeName = kwargs.pop('outputTreeName','AnalysisTree')
        self.shift = kwargs.pop('shift','')
        self.outputTreeName = outputTreeName
        if hasProgress:
            self.pbar = kwargs.pop('progressbar',ProgressBar(widgets=['{0}: '.format(outputTreeName),' ',SimpleProgress(),' events ',Percentage(),' ',Bar(),' ',ETA()]))
        # preselection
        if not hasattr(self,'preselection'): self.preselection = '1'
        # input files
        self.fileNames = []
        if isinstance(inputFileNames, basestring): # inputFiles is a file name
            if os.path.isfile(inputFileNames):     # single file
                if inputFileNames[-4:] == 'root':  # file is a root file
                    self.fileNames += [inputFileNames]
                else:                          # file is list of files
                    with open(inputFileNames,'r') as f:
                        for line in f:
                            self.fileNames += [line.strip()]
        else:
            self.fileNames = inputFileNames # already a python list or a cms.untracked.vstring()
        if not isinstance(outputFileName, basestring): # its a cms.string(), get value
            outputFileName = outputFileName.value()
        # input tchain
        self.treename = '{0}/{1}'.format(inputTreeDirectory,inputTreeName)
        tchainLumi = ROOT.TChain('{0}/{1}'.format(inputTreeDirectory,inputLumiName))
        self.totalEntries = 0
        logging.info('Getting Lumi information')
        #self.skims = {}
        for f,fName in enumerate(self.fileNames):
            if fName.startswith('/store'): fName = 'root://cmsxrootd.hep.wisc.edu//{0}'.format(fName)
            tfile = ROOT.TFile.Open(fName)
            tree = tfile.Get(self.treename)
            #skimName = 'skim{0}'.format(f)
            #tree.Draw('>>{0}'.format(skimName),self.preselection,'entrylist')
            #skimlist = ROOT.gDirectory.Get(skimName)
            #listEvents = skimlist.GetN()
            #self.skims[f] = skimlist
            #self.totalEntries += listEvents
            self.totalEntries += tree.GetEntries()
            if not hasattr(self,'version'):
                tree.GetEntry(1)
                if hasattr(tree,'provenance'):
                    ver = tree.provenance[0].split('_')
                    self.version = ''.join([ver[1],ver[2],'X'])
                else:
                    self.version = getCMSSWVersion()
            tfile.Close('R')
            tchainLumi.Add(fName)
        # get the lumi info
        self.numLumis = tchainLumi.GetEntries()
        self.numEvents = 0
        self.summedWeights = 0
        for entry in xrange(self.numLumis):
            tchainLumi.GetEntry(entry)
            self.numEvents += tchainLumi.nevents
            self.summedWeights += tchainLumi.summedWeights
        logging.info("Will process {0} lumi sections with {1} events ({2}).".format(self.numLumis,self.numEvents,self.summedWeights))
        self.flush()
        if not len(self.fileNames): raise Exception
        # other input files
        self.pileupWeights = PileupWeights(self.version)
        self.fakeRates = FakeRates(self.version)
        self.leptonScales = LeptonScales(self.version)
        self.triggerScales = TriggerScales(self.version)
        self.triggerPrescales = TriggerPrescales(self.version)
        self.zptGenWeight = ZptGenWeight(self.version)
        # tfile
        self.outfile = ROOT.TFile(outputFileName,"recreate")
        # cut tree
        self.cutTree = CutTree()
        # analysis tree
        self.tree = AnalysisTree(outputTreeName)
        self.eventsStored = 0

        # some things we always need:

        dysamples = [
            'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        ]


        # pileup
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[0], 'pileupWeight', 'F')
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[1], 'pileupWeightUp', 'F')
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[2], 'pileupWeightDown', 'F')
        self.tree.add(lambda cands: self.event.vertices_count(), 'numVertices', 'I')

        # gen
        self.tree.add(lambda cands: self.event.nTrueVertices(), 'numTrueVertices', 'I')
        self.tree.add(lambda cands: self.event.NUP(), 'NUP', 'I')
        self.tree.add(lambda cands: self.event.isData(), 'isData', 'I')
        if any([x in fName for x in dysamples]):
            self.tree.add(lambda cands: self.zptGenWeight.weight(self.gen), 'genWeight', 'F')
        else:
            self.tree.add(lambda cands: self.event.genWeight(), 'genWeight', 'F')
        self.tree.add(lambda cands: self.event.numGenJets(), 'numGenJets', 'I')
        self.tree.add(lambda cands: self.event.genHT(), 'genHT', 'I')
        # scale shifts
        weightMap = {
            0: {'muR':1.0, 'muF':1.0},
            1: {'muR':1.0, 'muF':2.0},
            2: {'muR':1.0, 'muF':0.5},
            3: {'muR':2.0, 'muF':1.0},
            4: {'muR':2.0, 'muF':2.0},
            5: {'muR':2.0, 'muF':0.5},
            6: {'muR':0.5, 'muF':1.0},
            7: {'muR':0.5, 'muF':2.0},
            8: {'muR':0.5, 'muF':0.5},
        }
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[0] if len(self.event.genWeights())>0 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[0]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[1] if len(self.event.genWeights())>1 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[1]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[2] if len(self.event.genWeights())>2 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[2]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[3] if len(self.event.genWeights())>3 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[3]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[4] if len(self.event.genWeights())>4 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[4]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[5] if len(self.event.genWeights())>5 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[5]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[6] if len(self.event.genWeights())>6 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[6]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[7] if len(self.event.genWeights())>7 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[7]), 'F')
        self.tree.add(lambda cands: 0. if self.event.isData() else self.event.genWeights()[8] if len(self.event.genWeights())>8 else 0., 'genWeight_muR{muR:3.1f}_muF{muF:3.1f}'.format(**weightMap[8]), 'F')


    def __exit__(self, type, value, traceback):
        self.finish()

    def __del__(self):
        self.finish()

    def finish(self):
        print ''
        logging.info('Finishing')
        logging.info('Writing {0} events'.format(self.eventsStored))
        self.outfile.cd()
        cutflowHist = ROOT.TH1F('summedWeights','summedWeights',1,0,1)
        cutflowHist.SetBinContent(1,self.summedWeights)
        self.outfile.Write()
        self.outfile.Close()
        self.leptonScales.finish()

    def flush(self):
        sys.stdout.flush()
        sys.stderr.flush()

    #############################
    ### primary analysis loop ###
    #############################
    def analyze(self):
        '''
        The primary analyzer loop.
        '''
        logging.info('Beginning Analysis')
        start = time.time()
        new = start
        old = start
        if hasProgress:
            self.pbar.maxval = self.totalEntries
            self.pbar.start()
            total = 0
            for f, fName in enumerate(self.fileNames):
                if fName.startswith('/store'): fName = 'root://cmsxrootd.hep.wisc.edu//{0}'.format(fName)
                tfile = ROOT.TFile.Open(fName,'READ')
                tree = tfile.Get(self.treename)
                #skimName = 'skim{0}'.format(f)
                #tree.Draw('>>{0}'.format(skimName),self.preselection,'entrylist')
                #skimlist = ROOT.gDirectory.Get(skimName)
                #listEvents = skimlist.GetN()
                #for r in xrange(listEvents):
                for row in tree:
                    total += 1
                    #tree.GetEntry(skimlist.Next())
                    self.pbar.update(total)
                    # load objects
                    self.event     = Event(tree)
                    if self.event.isData(): self.shift = ''
                    if not self.event.isData(): self.gen = [GenParticle(tree,entry=i) for i in range(tree.genParticles_count)]
                    self.electrons = [Electron(tree,entry=i,shift=self.shift) for i in range(tree.electrons_count)]
                    self.muons     = [Muon(tree,entry=i,shift=self.shift) for i in range(tree.muons_count)]
                    self.taus      = [Tau(tree,entry=i,shift=self.shift) for i in range(tree.taus_count)]
                    self.photons   = [Photon(tree,entry=i,shift=self.shift) for i in range(tree.photons_count)]
                    self.jets      = [Jet(tree,entry=i,shift=self.shift) for i in range(tree.jets_count)]
                    self.pfmet     = Met(tree,shift=self.shift)
                    # call per row action
                    self.perRowAction()
                tfile.Close('R')
            self.pbar.update(self.totalEntries)
        else:
            total = 0
            for f, fName in enumerate(self.fileNames):
                if fName.startswith('/store'): fName = 'root://cmsxrootd.hep.wisc.edu//{0}'.format(fName)
                logging.info('Processing file {0} of {1}'.format(f+1, len(self.fileNames)))
                tfile = ROOT.TFile.Open(fName,'READ')
                tree = tfile.Get(self.treename)
                #skimName = 'skim{0}'.format(f)
                #tree.Draw('>>{0}'.format(skimName),self.preselection,'entrylist')
                #skimlist = ROOT.gDirectory.Get(skimName)
                #listEvents = skimlist.GetN()
                #for r in xrange(listEvents):
                for row in tree:
                    total += 1
                    if total==2: start = time.time() # just ignore first event for timing
                    #tree.GetEntry(skimlist.Next())
                    if total % 1000 == 1:
                        cur = time.time()
                        elapsed = cur-start
                        remaining = float(elapsed)/total * float(self.totalEntries) - float(elapsed)
                        mins, secs = divmod(int(remaining),60)
                        hours, mins = divmod(mins,60)
                        logging.info('{0}: Processing event {1}/{2} - {3}:{4:02d}:{5:02d} remaining'.format(self.outputTreeName,total,self.totalEntries,hours,mins,secs))
                        self.flush()
                    # load objects
                    self.event     = Event(tree)
                    if self.event.isData(): self.shift = ''
                    if not self.event.isData(): self.gen = [GenParticle(tree,entry=i) for i in range(tree.genParticles_count)]
                    self.electrons = [Electron(tree,entry=i,shift=self.shift) for i in range(tree.electrons_count)]
                    self.muons     = [Muon(tree,entry=i,shift=self.shift) for i in range(tree.muons_count)]
                    self.taus      = [Tau(tree,entry=i,shift=self.shift) for i in range(tree.taus_count)]
                    self.photons   = [Photon(tree,entry=i,shift=self.shift) for i in range(tree.photons_count)]
                    self.jets      = [Jet(tree,entry=i,shift=self.shift) for i in range(tree.jets_count)]
                    self.pfmet     = Met(tree,shift=self.shift)
                    # call per row action
                    self.perRowAction()
                tfile.Close('R')

    def perRowAction(self):
        '''Per row action, can be overridden'''
        # select candidates
        cands = self.selectCandidates()
        cands['event'] = self.event

        # store event?
        goodToStore = self.cutTree.evaluate(cands)

        # do we store the tree?
        if not goodToStore: return

        self.tree.fill(cands)
        self.eventsStored += 1
        #self.outfile.Flush()

    def selectCandidates(self):
        '''
        Select candidates
            format should be:
            candidates = {
                "objectName" : ("collectionName", position),
                ...
            }
        '''
        logging.warning("You must override selectCandidates.")
        return {}

    #################
    ### utilities ###
    #################
    def findDecay(self,m_pdgid,d1_pdgid,d2_pdgid):
        '''Check if requested decay present in event'''
        for g in self.gen:
            if m_pdgid==g.pdgId():
                if (
                    (d1_pdgid==g.daughter_1()
                    and d2_pdgid==g.daughter_2())
                    or (d1_pdgid==g.daughter_2()
                    and d2_pdgid==g.daughter_1())
                   ):
                    return True
        return False

    def getCands(self,coll,func):
        cands = []
        for cand in coll:
            if func(cand): cands += [cand]
        return cands

    def cleanCands(self,src,other,dr):
        cleaned = []
        for s in src:
            keep = True
            for o in other:
                if deltaR(s.eta(),s.phi(),o.eta(),o.phi())<dr:
                    keep = False
            if keep:
                cleaned += [s]
        return cleaned

    def getCollectionString(self,cand):
        if isinstance(cand,Electron): return 'e'
        elif isinstance(cand,Muon):   return 'm'
        elif isinstance(cand,Tau):    return 't'
        elif isinstance(cand,Photon): return 'g'
        elif isinstance(cand,Jet):    return 'j'
        else:                         return 'a'

    ##########################
    ### add object to tree ###
    ##########################
    def addCandVar(self,label,varLabel,var,rootType):
        '''Add a variable for a cand'''
        self.tree.add(lambda cands: getattr(cands[label],var)(), '{0}_{1}'.format(label,varLabel), rootType)

    def addFlavorDependentCandVar(self,label,varLabel,varMap,rootType):
        '''Add a variable for a cand based on flavor'''
        self.tree.add(lambda cands: getattr(cands[label],varMap[cands[label].collName])() if cands[label].collName in varMap else 0., '{0}_{1}'.format(label,varLabel), rootType)

    def addMet(self,label):
        '''Add Met variables'''
        self.addCandVar(label,'pt','et','F')
        self.addCandVar(label,'phi','phi','F')

    def addJet(self,label):
        '''Add variables relevant for jets'''
        self.addCandVar(label,'pt','pt','F')
        self.addCandVar(label,'eta','eta','F')
        self.addCandVar(label,'phi','phi','F')
        self.addCandVar(label,'energy','energy','F')

    def addLepton(self,label):
        '''Add variables relevant for leptons'''
        self.addCandVar(label,'pt','pt','F')
        self.addCandVar(label,'eta','eta','F')
        self.addCandVar(label,'phi','phi','F')
        self.addCandVar(label,'energy','energy','F')
        self.addCandVar(label,'charge','charge','I')
        self.addCandVar(label,'dz','dz','F')
        self.addCandVar(label,'pdgId','pdgId','I')
        self.addFlavorDependentCandVar(label,'dxy',      {'electrons':'dB2D',          'muons':'dB2D', 'taus':'dxy',  '':''},'F')
        self.addFlavorDependentCandVar(label,'isolation',{'electrons':'relPFIsoRhoR03','muons':'relPFIsoDeltaBetaR04','':''},'F')
        self.addFlavorDependentCandVar(label,'genMatch',       {'electrons':'genMatch',       'muons':'genMatch',  'taus':'genJetMatch', '':''},'I')
        self.tree.add(lambda cands: self.genDeltaR(cands[label]) if isinstance(cands[label],Electron) or isinstance(cands[label],Muon) else self.genJetDeltaR(cands[label]), '{0}_genDeltaR'.format(label), 'F')
        self.addFlavorDependentCandVar(label,'genStatus',      {'electrons':'genStatus',      'muons':'genStatus', 'taus':'genJetStatus','':''},'I')
        self.addFlavorDependentCandVar(label,'genPdgId',       {'electrons':'genPdgId',       'muons':'genPdgId',  'taus':'genJetPdgId', '':''},'I')
        self.addFlavorDependentCandVar(label,'genPt',          {'electrons':'genPt',          'muons':'genPt',     'taus':'genJetPt',    '':''},'F')
        self.addFlavorDependentCandVar(label,'genEta',         {'electrons':'genEta',         'muons':'genEta',    'taus':'genJetEta',   '':''},'F')
        self.addFlavorDependentCandVar(label,'genPhi',         {'electrons':'genPhi',         'muons':'genPhi',    'taus':'genJetPhi',   '':''},'F')
        self.addFlavorDependentCandVar(label,'genEnergy',      {'electrons':'genEnergy',      'muons':'genEnergy', 'taus':'genJetEnergy','':''},'F')
        self.addFlavorDependentCandVar(label,'genCharge',      {'electrons':'genCharge',      'muons':'genCharge', 'taus':'genJetCharge','':''},'I')
        self.addFlavorDependentCandVar(label,'genIsPrompt',    {'electrons':'genIsPrompt',    'muons':'genIsPrompt',    '':''},'I')
        self.addFlavorDependentCandVar(label,'genIsFromTau',   {'electrons':'genIsFromTau',   'muons':'genIsFromTau',   '':''},'I')
        self.addFlavorDependentCandVar(label,'genIsFromHadron',{'electrons':'genIsFromHadron','muons':'genIsFromHadron','':''},'I')

    def genDeltaR(self,cand):
        '''Get the gen level deltaR'''
        if cand.genMatch()==0: return 0.
        eta = cand.eta()
        genEta = cand.genEta()
        phi = cand.phi()
        genPhi = cand.genPhi()
        return deltaR(eta,phi,genEta,genPhi)

    def genJetDeltaR(self,cand):
        '''Get the gen level deltaR'''
        if cand.genJetMatch()==0: return 0.
        eta = cand.eta()
        genEta = cand.genJetEta()
        phi = cand.phi()
        genPhi = cand.genJetPhi()
        return deltaR(eta,phi,genEta,genPhi)

    def addDiJet(self,label):
        '''Add variables relevant for a dijet candidate'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaR','deltaR','F')
        self.addCandVar(label,'deltaEta','deltaEta','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')
        self.addCandVar(label,'energy','Energy','F')

    def addDiLepton(self,label):
        '''Add variables relevant for a dilepton candidate'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaR','deltaR','F')
        self.addCandVar(label,'deltaEta','deltaEta','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')
        self.addCandVar(label,'energy','Energy','F')

    def addLeptonMet(self,label):
        '''Add variables related to a lepton + met'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')

    def addComposite(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'energy','Energy','F')

    def addCompositeMet(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')

