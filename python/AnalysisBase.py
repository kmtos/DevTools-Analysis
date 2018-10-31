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
from ZZGenWeight import ZZGenWeight

from utilities import deltaR, deltaPhi
from DevTools.Utilities.utilities import getCMSSWVersion
from Candidates import *
from leptonId import passHppLoose, passHppMedium, passHppTight
from photonId import passPhoton, passPreselection, passPreselectionNoElectronVeto, passElectronVeto

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
        self.events = []
        self.outputTreeName = outputTreeName
        if hasProgress:
            self.pbar = kwargs.pop('progressbar',ProgressBar(widgets=['{0}: '.format(outputTreeName),' ',SimpleProgress(),' events ',Percentage(),' ',Bar(),' ',ETA()]))
        # preselection
        if not hasattr(self,'preselection'): self.preselection = '1'
        # input files
        self.fileNames = []
        if os.path.isfile('PSet.py'):                # grab input files from crab pset
            import PSet
            self.fileNames = list(PSet.process.source.fileNames)
        elif isinstance(inputFileNames, basestring): # inputFiles is a file name
            if os.path.isfile(inputFileNames):       # single file
                if inputFileNames[-4:] == 'root':    # file is a root file
                    self.fileNames += [inputFileNames]
                else:                                # file is list of files
                    with open(inputFileNames,'r') as f:
                        for line in f:
                            self.fileNames += [line.strip()]
        else:
            self.fileNames = inputFileNames          # already a python list or a cms.untracked.vstring()
        if not isinstance(outputFileName, basestring): # its a cms.string(), get value
            outputFileName = outputFileName.value()
        # test for hdfs
        #self.hasHDFS = os.path.exists('/hdfs/store/user')
        self.hasHDFS = False
        # input tchain
        self.treename = '{0}/{1}'.format(inputTreeDirectory,inputTreeName)
        luminame = '{0}/{1}'.format(inputTreeDirectory,inputLumiName) if inputLumiName else ''
        #tchainLumi = ROOT.TChain(luminame)
        self.totalEntries = 0
        self.numLumis = 0
        self.numEvents = 0
        self.summedWeights = 0
        logging.info('Getting Lumi information')
        #self.skims = {}
        for f,fName in enumerate(self.fileNames):
            if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
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
            if luminame:
                lumitree = tfile.Get(luminame)
                for entry in lumitree:
                    self.numLumis += 1
                    self.numEvents += lumitree.nevents
                    self.summedWeights += lumitree.summedWeights
            tfile.Close('R')
            #tchainLumi.Add(fName)
        # get the lumi info
        #self.numLumis = tchainLumi.GetEntries()
        #self.numEvents = 0
        #self.summedWeights = 0
        #for entry in xrange(self.numLumis):
        #    tchainLumi.GetEntry(entry)
        #    self.numEvents += tchainLumi.nevents
        #    self.summedWeights += tchainLumi.summedWeights
        logging.info('Analysis is running with version {0}'.format(self.version))
        logging.info("Will process {0} lumi sections with {1} events ({2}).".format(self.numLumis,self.numEvents,self.summedWeights))
        self.flush()
        if not len(self.fileNames): raise Exception
        # other input files
        if 'UpsilonMuMu_UpsilonPt6_TuneCUEP8M1_13TeV-pythia8-evtgen' in fName:
            self.pileupWeights = PileupWeights(self.version,profile='PUSpring16')
        else:
            self.pileupWeights = PileupWeights(self.version)
        self.fakeRates = FakeRates(self.version)
        self.leptonScales = LeptonScales(self.version)
        self.triggerScales = TriggerScales(self.version)
        self.triggerPrescales = TriggerPrescales(self.version)
        self.zptGenWeight = ZptGenWeight(self.version)
        self.zzGenWeight = ZZGenWeight(self.version)
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

        qqzzsamples = [
            'ZZTo4L_13TeV_powheg_pythia8',
            'ZZTo4L_13TeV-amcatnloFXFX-pythia8',
        ]

        # pileup
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[0], 'pileupWeight', 'F')
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[1], 'pileupWeightUp', 'F')
        self.tree.add(lambda cands: self.pileupWeights.weight(self.event)[2], 'pileupWeightDown', 'F')
        self.tree.add(lambda cands: self.event.vertices_count(), 'numVertices', 'I')
        self.tree.add(lambda cands: self.event.rho(), 'rho', 'F')

        # gen
        self.tree.add(lambda cands: self.event.nTrueVertices(), 'numTrueVertices', 'I')
        self.tree.add(lambda cands: self.event.NUP(), 'NUP', 'I')
        self.tree.add(lambda cands: self.event.isData(), 'isData', 'I')
        self.tree.add(lambda cands: self.event.genWeight(), 'genWeight', 'F')
        if any([x in fName for x in dysamples]):
            self.tree.add(lambda cands: self.zptGenWeight.weight(self.gen), 'zPtWeight', 'F')
        if any([x in fName for x in qqzzsamples]):
            self.tree.add(lambda cands: self.zzGenWeight.weight(self.gen), 'qqZZkfactor', 'F')
        self.tree.add(lambda cands: self.event.numGenJets(), 'numGenJets', 'I')
        self.tree.add(lambda cands: self.event.genHT(), 'genHT', 'I')
        #self.tree.add(lambda cands: 0 if self.event.isData() else self.event.ptHat() if hasattr(self.event,'ptHat') else 0, 'ptHat', 'F')
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
        cutflowHist = ROOT.TH1D('summedWeights','summedWeights',1,0,1)
        cutflowHist.SetBinContent(1,self.summedWeights)
        self.outfile.Write()
        self.outfile.Close()
        self.leptonScales.finish()

    def flush(self):
        sys.stdout.flush()
        sys.stderr.flush()

    def get_event(self):
        return '{run}:{lumi}:{event}'.format(run=self.event.run(),lumi=self.event.lumi(),event=self.event.event())

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
                if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
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
                    self.setupEvent(tree)
                    self.perRowAction()
                tfile.Close('R')
            self.pbar.update(self.totalEntries)
        else:
            total = 0
            for f, fName in enumerate(self.fileNames):
                if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
                logging.info('Processing file {0} of {1}: {2}'.format(f+1, len(self.fileNames), fName))
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
                    self.setupEvent(tree)
                    self.perRowAction()
                tfile.Close('R')

    def setupEvent(self,tree):
        '''Setup the event objects'''
        # load objects
        self.event     = Event(tree)
        if self.event.isData(): self.shift = ''
        if not self.event.isData(): self.gen = [GenParticle(tree,entry=i) for i in range(tree.genParticles_count)]
        self.electrons = [Electron(tree,entry=i,shift=self.shift) for i in range(tree.electrons_count)]
        self.muons     = [Muon(tree,entry=i,shift=self.shift) for i in range(tree.muons_count)]
        self.taus      = [Tau(tree,entry=i,shift=self.shift) for i in range(tree.taus_count)]
        if hasattr(tree, 'photons_count'): self.photons   = [Photon(tree,entry=i,shift=self.shift) for i in range(tree.photons_count)]
        self.jets      = [Jet(tree,entry=i,shift=self.shift) for i in range(tree.jets_count)]
        self.pfmet     = Met(tree,shift=self.shift)
        if hasattr(tree, 'packed_count'): self.packed   = [Candidate(tree,entry=i,collName='packed') for i in range(tree.packed_count)]

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
                    break
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

    def checkTrigger(self,*datasets,**triggerNames):
        '''Check trigger using trigger map'''
        isData = self.event.isData()>0.5
        # reject triggers if they are in another dataset
        # looks for the dataset name in the filename
        # for MC it accepts any
        reject = True if isData else False
        for dataset in datasets:
            # if we match to the dataset, start accepting triggers
            if dataset in self.fileNames[0] and isData: reject = False
            for trigger in triggerNames[dataset]:
                var = '{0}Pass'.format(trigger)
                passTrigger = getattr(self.event,var)()
                if passTrigger>0.5:
                    # it passed the trigger
                    # in data: reject if it corresponds to a higher dataset
                    return False if reject else True
            # dont check the rest of data
            if dataset in self.fileNames[0] and isData: break
        return False

    def metFilter(self,cands):
        if not self.event.isData(): return True
        filterList = [
            'HBHENoiseFilter',
            'HBHENoiseIsoFilter',
            'globalTightHalo2016Filter',
            'EcalDeadCellTriggerPrimitiveFilter',
            'goodVertices',
            'eeBadScFilter',
            'noBadMuons',
            'BadChargedCandidateFilter',
        ]
        notFilterList = [
            # Dont use, "noBadMuons" contains this
            #'duplicateMuons',
            #'badMuons',
        ]
        for f in filterList:
            if getattr(self.event,f)()==0:
                logging.info('Rejecting event {0}:{1}:{2} for {3}={4}'.format(self.event.run(), self.event.lumi(), self.event.event(), f, getattr(self.event,f)()))
                self.report_failure('fails {}'.format(f))
                return False
        for f in notFilterList:
            if getattr(self.event,f)()>0:
                logging.info('Rejecting event {0}:{1}:{2} for {3}={4}'.format(self.event.run(), self.event.lumi(), self.event.event(), f, getattr(self.event,f)()))
                self.report_failure('fails {}'.format(f))
                return False
        return True

    def report_failure(self,message):
        if self.events:
            event = self.get_event()
            if event in self.events:
                 print event, message

    ##################
    ### Common IDs ###
    ##################
    # override in derived class if desired
    # base selections for H++ analysis
    def passLooseNew(self,cand):
        return passHppLoose(cand,new=True)

    def passLoose(self,cand):
        return passHppLoose(cand)

    def passMedium(self,cand):
        return passHppMedium(cand)

    def passTight(self,cand):
        return passHppTight(cand)

    def passPhotonId(self,cand):
        return passPhoton(cand)

    def passElectronVeto(self,cand):
        return passElectronVeto(cand)

    def passPhotonPreselection(self,cand):
        return passPreselection(cand)

    def passPhotonPreselectionNoElectronVeto(self,cand):
        return passPreselectionNoElectronVeto(cand)

    def looseScale(self,cand):
        #key = 'CutbasedVeto' if abs(cand.eta())<1.479 else 'CutbasedLoose'
        if cand.collName=='muons': return self.leptonScales.getScale('MediumIDLooseIso',cand,doError=True)
        #elif cand.collName=='electrons': return self.leptonScales.getScale(key,cand,doError=True)
        elif cand.collName=='electrons': return self.leptonScales.getScale('CutbasedLoose',cand,doError=True)
        elif cand.collName=='taus': return self.leptonScales.getScale('',cand,doError=True)
        else: return [1.,1.,1.]

    def mediumScale(self,cand):
        if cand.collName=='muons': return self.leptonScales.getScale('MediumIDTightIso',cand,doError=True)
        elif cand.collName=='electrons': return self.leptonScales.getScale('CutbasedMedium',cand,doError=True)
        elif cand.collName=='taus': return self.leptonScales.getScale('',cand,doError=True)
        else: return [1.,1.,1.]

    def tightScale(self,cand):
        #if cand.collName=='muons': return self.leptonScales.getScale('MediumIDTightIso',cand,doError=True)
        if cand.collName=='muons': return self.leptonScales.getScale('TightIDTightIso',cand,doError=True)
        elif cand.collName=='electrons': return self.leptonScales.getScale('CutbasedTight',cand,doError=True)
        elif cand.collName=='taus': return self.leptonScales.getScale('',cand,doError=True)
        else: return [1.,1.,1.]

    def photonPreselectionScale(self,cand):
        return self.leptonScales.getScale('Preselection',cand,doError=True)

    def photonMVAScale(self,cand):
        return self.leptonScales.getScale('MVA0p0',cand,doError=True)

    def photonMVAPreselectionScale(self,cand):
        return self.leptonScales.getScale('MVA0p0Pre',cand,doError=True)

    def mediumNewDMFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppMedium','HppLooseNew',doError=True,doDM=True)

    def tightNewDMFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppLooseNew',doError=True,doDM=True)

    def mediumDMFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppMedium','HppLoose',doError=True,doDM=True)

    def tightDMFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppLoose',doError=True,doDM=True)

    def mediumNewFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppMedium','HppLooseNew',doError=True)

    def tightNewFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppLooseNew',doError=True)

    def mediumFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppMedium','HppLoose',doError=True)

    def tightFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppLoose',doError=True)

    def tightFromMediumFakeRate(self,cand):
        return self.fakeRates.getFakeRate(cand,'HppTight','HppMedium',doError=True)

    def getPassingCands(self,mode,*colls):
        if mode=='Loose': passMode = self.passLoose
        elif mode=='LooseNew': passMode = self.passLooseNew
        elif mode=='Medium': passMode = self.passMedium
        elif mode=='Tight': passMode = self.passTight
        elif mode=='Photon': passMode = self.passPhotonId
        elif mode=='ElectronVeto': passMode = self.passElectronVeto
        elif mode=='PhotonPreselection': passMode = self.passPhotonPreselection
        elif mode=='PhotonPreselectionNoElectronVeto': passMode = self.passPhotonPreselectionNoElectronVeto
        else: return []
        cands = []
        for coll in colls:
            cands += self.getCands(coll,passMode)
        return cands

    def numJets(self,cleanJets,mode,pt):
        jetColl = self.getCands(
            cleanJets,
            lambda cand: getattr(cand,mode)()>0.5 and cand.pt()>pt
        )
        return len(jetColl)

    ##########################
    ### add object to tree ###
    ##########################
    def addTriggers(self):
        if self.version=='76X':
            # single lepton
            self.tree.add(lambda cands: self.event.IsoMu20Pass(), 'pass_IsoMu20', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu20Pass(), 'pass_IsoTkMu20', 'I')
            self.tree.add(lambda cands: self.event.Ele23_WPLoose_GsfPass(), 'pass_Ele23_WPLoose_Gsf', 'I')
            # double lepton
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.DoubleMediumIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg', 'I')
        else:
            # single lepton
            self.tree.add(lambda cands: self.event.IsoMu24Pass(), 'pass_IsoMu24', 'I')
            self.tree.add(lambda cands: self.event.IsoTkMu24Pass(), 'pass_IsoTkMu24', 'I')
            self.tree.add(lambda cands: self.event.Mu50Pass(), 'pass_Mu50', 'I')
            self.tree.add(lambda cands: self.event.TkMu50Pass(), 'pass_TkMu50', 'I')
            self.tree.add(lambda cands: self.event.Ele27_WPTight_GsfPass(), 'pass_Ele27_WPTight_Gsf', 'I')
            # double lepton
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZPass(), 'pass_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu17_Mu8_SameSign_DZPass(), 'pass_Mu17_Mu8_SameSign_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu20_Mu10_SameSign_DZPass(), 'pass_Mu20_Mu10_SameSign_DZ', 'I')
            self.tree.add(lambda cands: self.event.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVLPass(), 'pass_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL', 'I')
            self.tree.add(lambda cands: self.event.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZPass(), 'pass_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ', 'I')
            self.tree.add(lambda cands: self.event.DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_RegPass(), 'pass_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg', 'I')
            # triple lepton
            self.tree.add(lambda cands: self.event.Ele16_Ele12_Ele8_CaloIdL_TrackIdLPass(), 'pass_Ele16_Ele12_Ele8_CaloIdL_TrackIdL', 'I')
            self.tree.add(lambda cands: self.event.Mu8_DiEle12_CaloIdL_TrackIdLPass(), 'pass_Mu8_DiEle12_CaloIdL_TrackIdL', 'I')
            self.tree.add(lambda cands: self.event.DiMu9_Ele9_CaloIdL_TrackIdLPass(), 'pass_DiMu9_Ele9_CaloIdL_TrackIdL', 'I')
            self.tree.add(lambda cands: self.event.TripleMu_12_10_5Pass(), 'pass_TripleMu_12_10_5', 'I')

    def addPhotonTriggers(self):
        self.tree.add(lambda cands: self.event.DoublePhoton60Pass(), 'pass_DoublePhoton60', 'I')
        self.tree.add(lambda cands: self.event.Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90Pass(), 'pass_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90', 'I')
        self.tree.add(lambda cands: self.event.Photon175Pass(), 'pass_Photon175', 'I')
        

    def addCandVar(self,label,varLabel,var,rootType,permissive=False):
        '''Add a variable for a cand'''
        if permissive:
            self.tree.add(lambda cands: getattr(cands[label],var)() if cands[label].has(var) else 0., '{0}_{1}'.format(label,varLabel), rootType)
        else:
            self.tree.add(lambda cands: getattr(cands[label],var)(), '{0}_{1}'.format(label,varLabel), rootType)

    def addFlavorDependentCandVar(self,label,varLabel,varMap,rootType,default=0):
        '''Add a variable for a cand based on flavor'''
        self.tree.add(lambda cands: getattr(cands[label],varMap[cands[label].collName])() if cands[label].collName in varMap else default, '{0}_{1}'.format(label,varLabel), rootType)

    def addMet(self,label):
        '''Add Met variables'''
        self.addCandVar(label,'pt','et','F')
        self.addCandVar(label,'phi','phi','F')
        self.addCandVar(label,'cov00','cov00','F',permissive=True)
        self.addCandVar(label,'cov01','cov01','F',permissive=True)
        self.addCandVar(label,'cov10','cov10','F',permissive=True)
        self.addCandVar(label,'cov11','cov11','F',permissive=True)

    def addCandidate(self,label):
        '''Add variables relevant for all objects'''
        self.addCandVar(label,'pt','pt','F')
        self.addCandVar(label,'eta','eta','F')
        self.addCandVar(label,'phi','phi','F')
        self.addCandVar(label,'energy','energy','F')
        self.addCandVar(label,'mass','mass','F')

    def addJet(self,label):
        '''Add variables relevant for jets'''
        self.addCandidate(label)
        self.addCandVar(label,'CSVv2','pfCombinedInclusiveSecondaryVertexV2BJetTags','F')
        self.addCandVar(label,'hadronFlavour','hadronFlavour','I')
        self.addCandVar(label,'passCSVv2T','passCSVv2T','I')
        self.addCandVar(label,'passCSVv2M','passCSVv2M','I')
        self.addCandVar(label,'passCSVv2L','passCSVv2L','I')

    def addLepton(self,label,doId=False,doScales=False,doFakes=False,doErrors=False):
        '''Add variables relevant for leptons'''
        self.addCandidate(label)
        self.addCandVar(label,'charge','charge','I')
        self.addCandVar(label,'dz','dz','F')
        self.addCandVar(label,'pdgId','pdgId','I')
        self.addFlavorDependentCandVar(label,'dxy',            {'electrons':'dB2D',           'muons':'dB2D',      'taus':'dxy',         '':''},'F')
        self.addFlavorDependentCandVar(label,'isolation',      {'electrons':'relPFIsoRhoR03', 'muons':'relPFIsoDeltaBetaR04',            '':''},'F')
        self.addCandVar(label,'genMatch','genMatch','I')
        self.tree.add(lambda cands: self.genDeltaR(cands[label]), '{0}_genDeltaR'.format(label), 'F')
        self.addCandVar(label,'genStatus','genStatus','F')
        self.addCandVar(label,'genPdgId','genPdgId','F')
        self.addCandVar(label,'genPt','genPt','F')
        self.addCandVar(label,'genEta','genEta','F')
        self.addCandVar(label,'genPhi','genPhi','F')
        self.addCandVar(label,'genEnergy','genEnergy','F')
        self.addCandVar(label,'genCharge','genCharge','F')
        self.addFlavorDependentCandVar(label,'decayMode',         {'taus':'decayMode', '':''},                         'I', default=-1)
        self.addFlavorDependentCandVar(label,'decayModeFinding',  {'taus':'decayModeFinding', '':''},                  'I', default=-1)
        self.addFlavorDependentCandVar(label,'tauMVAold',         {'taus':'byIsolationMVArun2v1DBoldDMwLTraw', '':''}, 'F', default=-2)
        self.addFlavorDependentCandVar(label,'tauMVAnew',         {'taus':'byIsolationMVArun2v1DBnewDMwLTraw', '':''}, 'F', default=-2)
        self.addFlavorDependentCandVar(label,'genJetMatch',       {'taus':'genJetMatch', '':''},                       'I', default=-1)
        self.tree.add(lambda cands: 0 if isinstance(cands[label],Electron) or isinstance(cands[label],Muon) else self.genJetDeltaR(cands[label]), '{0}_genJetDeltaR'.format(label), 'F')
        self.addFlavorDependentCandVar(label,'genJetStatus',      {'taus':'genJetStatus','':''},'I')
        self.addFlavorDependentCandVar(label,'genJetPdgId',       {'taus':'genJetPdgId', '':''},'I')
        self.addFlavorDependentCandVar(label,'genJetPt',          {'taus':'genJetPt',    '':''},'F')
        self.addFlavorDependentCandVar(label,'genJetEta',         {'taus':'genJetEta',   '':''},'F')
        self.addFlavorDependentCandVar(label,'genJetPhi',         {'taus':'genJetPhi',   '':''},'F')
        self.addFlavorDependentCandVar(label,'genJetEnergy',      {'taus':'genJetEnergy','':''},'F')
        self.addFlavorDependentCandVar(label,'genJetCharge',      {'taus':'genJetCharge','':''},'I')
        self.addFlavorDependentCandVar(label,'genIsPrompt',       {'electrons':'genIsPrompt',    'muons':'genIsPrompt',      '':''}, 'I', default=-1)
        self.addFlavorDependentCandVar(label,'genIsFromTau',      {'electrons':'genIsFromTau',   'muons':'genIsFromTau',     '':''}, 'I', default=-1)
        self.addFlavorDependentCandVar(label,'genIsFromHadron',   {'electrons':'genIsFromHadron','muons':'genIsFromHadron',  '':''}, 'I', default=-1)
        if doId:
            self.tree.add(lambda cands: self.passLoose(cands[label]),                      '{0}_passLoose'.format(label), 'I')
            self.tree.add(lambda cands: self.passLooseNew(cands[label]),                   '{0}_passLooseNew'.format(label), 'I')
            self.tree.add(lambda cands: self.passMedium(cands[label]),                     '{0}_passMedium'.format(label), 'I')
            self.tree.add(lambda cands: self.passTight(cands[label]),                      '{0}_passTight'.format(label), 'I')
        if doScales:
            self.tree.add(lambda cands: self.looseScale(cands[label])[0],                  '{0}_looseScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.looseScale(cands[label])[1],     '{0}_looseScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.looseScale(cands[label])[2],     '{0}_looseScaleDown'.format(label), 'F')
            self.tree.add(lambda cands: self.mediumScale(cands[label])[0],                 '{0}_mediumScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumScale(cands[label])[1],    '{0}_mediumScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumScale(cands[label])[2],    '{0}_mediumScaleDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightScale(cands[label])[0],                  '{0}_tightScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightScale(cands[label])[1],     '{0}_tightScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightScale(cands[label])[2],     '{0}_tightScaleDown'.format(label), 'F')
        if doFakes:
            self.tree.add(lambda cands: self.mediumNewFakeRate(cands[label])[0],                     '{0}_mediumNewFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumNewFakeRate(cands[label])[1],        '{0}_mediumNewFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumNewFakeRate(cands[label])[2],        '{0}_mediumNewFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightNewFakeRate(cands[label])[0],                      '{0}_tightNewFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightNewFakeRate(cands[label])[1],         '{0}_tightNewFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightNewFakeRate(cands[label])[2],         '{0}_tightNewFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.mediumNewDMFakeRate(cands[label])[0],                   '{0}_mediumNewDMFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumNewDMFakeRate(cands[label])[1],      '{0}_mediumNewDMFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumNewDMFakeRate(cands[label])[2],      '{0}_mediumNewDMFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightNewDMFakeRate(cands[label])[0],                    '{0}_tightNewDMFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightNewDMFakeRate(cands[label])[1],       '{0}_tightNewDMFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightNewDMFakeRate(cands[label])[2],       '{0}_tightNewDMFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.mediumDMFakeRate(cands[label])[0],                      '{0}_mediumDMFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumDMFakeRate(cands[label])[1],         '{0}_mediumDMFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumDMFakeRate(cands[label])[2],         '{0}_mediumDMFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightDMFakeRate(cands[label])[0],                       '{0}_tightDMFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightDMFakeRate(cands[label])[1],          '{0}_tightDMFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightDMFakeRate(cands[label])[2],          '{0}_tightDMFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.mediumFakeRate(cands[label])[0],                        '{0}_mediumFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumFakeRate(cands[label])[1],           '{0}_mediumFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.mediumFakeRate(cands[label])[2],           '{0}_mediumFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightFakeRate(cands[label])[0],                         '{0}_tightFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightFakeRate(cands[label])[1],            '{0}_tightFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightFakeRate(cands[label])[2],            '{0}_tightFakeRateDown'.format(label), 'F')
            self.tree.add(lambda cands: self.tightFromMediumFakeRate(cands[label])[0],               '{0}_tightFromMediumFakeRate'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightFromMediumFakeRate(cands[label])[1],  '{0}_tightFromMediumFakeRateUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.tightFromMediumFakeRate(cands[label])[2],  '{0}_tightFromMediumFakeRateDown'.format(label), 'F')

    def addDetailedMuon(self,label):
        '''Add detailed  variables'''
        self.addCandVar(label,'isLooseMuon','isLooseMuon','I')
        self.addCandVar(label,'isMediumMuon','isMediumMuon','I')
        self.addCandVar(label,'isMediumMuonICHEP','isMediumMuonICHEP','I')
        self.addCandVar(label,'isTightMuon','isTightMuon','I')
        self.addCandVar(label,'isPFMuon','isPFMuon','I')
        self.addCandVar(label,'isGlobalMuon','isGlobalMuon','I')
        self.addCandVar(label,'isTrackerMuon','isTrackerMuon','I')
        self.addCandVar(label,'muonBestTrackType','muonBestTrackType','I')
        self.addCandVar(label,'segmentCompatibility','segmentCompatibility','F')
        self.addCandVar(label,'isGoodMuon','isGoodMuon','I')
        self.addCandVar(label,'highPurityTrack','highPurityTrack','I')
        self.addCandVar(label,'matchedStations','matchedStations','I')
        self.addCandVar(label,'validMuonHits','validMuonHits','I')
        self.addCandVar(label,'normalizedChi2','normalizedChi2','F')
        self.addCandVar(label,'validPixelHits','validPixelHits','I')
        self.addCandVar(label,'trackerLayers','trackerLayers','I')
        self.addCandVar(label,'pixelLayers','pixelLayers','I')
        self.addCandVar(label,'validTrackerFraction','validTrackerFraction','F')
        self.addCandVar(label,'bestTrackPtError','bestTrackPtError','F')
        self.addCandVar(label,'bestTrackPt','bestTrackPt','F')
        self.addCandVar(label,'trackerStandaloneMatch','trackerStandaloneMatch','F')
        self.addCandVar(label,'relPFIsoDeltaBetaR04','relPFIsoDeltaBetaR04','F')
        self.tree.add(lambda cands: cands[label].trackIso()/cands[label].pt() if cands[label].pt() else 0., '{0}_trackRelIso'.format(label), 'F')

    def addDetailedTau(self,label):
        '''Add detailed variables'''
        # against muon
        self.addCandVar(label,'againstMuonLoose3','againstMuonLoose3','I')
        self.addCandVar(label,'againstMuonTight3','againstMuonTight3','I')
        # against electron
        self.addCandVar(label,'againstElectronVLooseMVA6','againstElectronVLooseMVA6','I')
        self.addCandVar(label,'againstElectronLooseMVA6','againstElectronLooseMVA6','I')
        self.addCandVar(label,'againstElectronMediumMVA6','againstElectronMediumMVA6','I')
        self.addCandVar(label,'againstElectronTightMVA6','againstElectronTightMVA6','I')
        self.addCandVar(label,'againstElectronVTightMVA6','againstElectronVTightMVA6','I')
        #self.addCandVar(label,'decayModeFinding','decayModeFinding','I')
        #self.addCandVar(label,'decayMode','decayMode','I')
        # MVA
        self.addCandVar(label,'byIsolationMVArun2v1DBoldDMwLTraw','byIsolationMVArun2v1DBoldDMwLTraw','F')
        self.addCandVar(label,'byVLooseIsolationMVArun2v1DBoldDMwLT','byVLooseIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byLooseIsolationMVArun2v1DBoldDMwLT','byLooseIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byMediumIsolationMVArun2v1DBoldDMwLT','byMediumIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byTightIsolationMVArun2v1DBoldDMwLT','byTightIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byVTightIsolationMVArun2v1DBoldDMwLT','byVTightIsolationMVArun2v1DBoldDMwLT','I')
        self.addCandVar(label,'byIsolationMVArun2v1DBnewDMwLTraw','byIsolationMVArun2v1DBnewDMwLTraw','F')
        self.addCandVar(label,'byVLooseIsolationMVArun2v1DBnewDMwLT','byVLooseIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byLooseIsolationMVArun2v1DBnewDMwLT','byLooseIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byMediumIsolationMVArun2v1DBnewDMwLT','byMediumIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byTightIsolationMVArun2v1DBnewDMwLT','byTightIsolationMVArun2v1DBnewDMwLT','I')
        self.addCandVar(label,'byVTightIsolationMVArun2v1DBnewDMwLT','byVTightIsolationMVArun2v1DBnewDMwLT','I')
        # cutbased
        self.addCandVar(label,'byCombinedIsolationDeltaBetaCorrRaw3Hits','byCombinedIsolationDeltaBetaCorrRaw3Hits','F')
        self.addCandVar(label,'chargedIsoPtSum','chargedIsoPtSum','F')
        self.addCandVar(label,'neutralIsoPtSum','neutralIsoPtSum','F')
        self.addCandVar(label,'puCorrPtSum','puCorrPtSum','F')
        self.addCandVar(label,'footprintCorrection','footprintCorrection','F')
        self.addCandVar(label,'photonPtSumOutsideSignalCone','photonPtSumOutsideSignalCone','F')

    def addPhoton(self,label,doId=False,doScales=False,doFakes=False,doErrors=False):
        '''Add variables relevant for photons'''
        self.addCandidate(label)
        self.addCandVar(label,'mvaNonTrigValues','mvaNonTrigValues','F')
        self.addCandVar(label,'mvaNonTrigCategories','mvaNonTrigCategories','F')
        self.addCandVar(label,'r9','r9','F')
        self.addCandVar(label,'chargedIso','phoChargedIsolation','F')
        self.addCandVar(label,'neutralHadronIso','phoNeutralHadronIsolation','F')
        self.addCandVar(label,'photonIso','phoPhotonIsolation','F')
        if doId:
            self.tree.add(lambda cands: self.passPhotonId(cands[label]),           '{0}_passId'.format(label), 'I')
            self.tree.add(lambda cands: self.passElectronVeto(cands[label]),       '{0}_passElectronVeto'.format(label), 'I')
            self.tree.add(lambda cands: self.passPhotonPreselection(cands[label]), '{0}_passPreselection'.format(label), 'I')
            self.tree.add(lambda cands: self.passPhotonPreselectionNoElectronVeto(cands[label]), '{0}_passPreselectionNoElectronVeto'.format(label), 'I')
        if doScales:
            self.tree.add(lambda cands: self.photonPreselectionScale(cands[label])[0],                  '{0}_preselectionScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonPreselectionScale(cands[label])[1],     '{0}_preselectionScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonPreselectionScale(cands[label])[2],     '{0}_preselectionScaleDown'.format(label), 'F')
            self.tree.add(lambda cands: self.photonMVAScale(cands[label])[0],                           '{0}_mvaScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonMVAScale(cands[label])[1],              '{0}_mvaScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonMVAScale(cands[label])[2],              '{0}_mvaScaleDown'.format(label), 'F')
            self.tree.add(lambda cands: self.photonMVAPreselectionScale(cands[label])[0],               '{0}_mvaPreselectionScale'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonMVAPreselectionScale(cands[label])[1],  '{0}_mvaPreselectionScaleUp'.format(label), 'F')
            if doErrors: self.tree.add(lambda cands: self.photonMVAPreselectionScale(cands[label])[2],  '{0}_mvaPreselectionScaleDown'.format(label), 'F')

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

    def addDiCandidate(self,label):
        '''Add variables relevant for a two object candidate'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'M2','M2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaR','deltaR','F')
        self.addCandVar(label,'deltaEta','deltaEta','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')
        self.addCandVar(label,'energy','Energy','F')

    def addDiJet(self,label):
        '''Add variables relevant for a dijet candidate'''
        self.addDiCandidate(label)

    def addDiLepton(self,label):
        '''Add variables relevant for a dilepton candidate'''
        self.addDiCandidate(label)

    def addLeptonMet(self,label):
        '''Add variables related to a lepton + met'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'mt2','Mt2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'deltaPhi','deltaPhi','F')

    def addComposite(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mass','M','F')
        self.addCandVar(label,'M2','M2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
        self.addCandVar(label,'energy','Energy','F')

    def addCompositeMet(self,label):
        '''Add variables related to multi object variables'''
        self.addCandVar(label,'mt','Mt','F')
        self.addCandVar(label,'mt2','Mt2','F')
        self.addCandVar(label,'pt','Pt','F')
        self.addCandVar(label,'eta','Eta','F')
        self.addCandVar(label,'phi','Phi','F')
