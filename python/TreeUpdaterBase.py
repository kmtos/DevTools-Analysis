import logging
import os
import sys
import math
import time

import ROOT
ROOT.gROOT.SetBatch(True)

from PileupWeights import PileupWeights
from FakeRates import FakeRates
from LeptonScales import LeptonScales
from TriggerScales import TriggerScales
from TriggerPrescales import TriggerPrescales
from ZptGenWeight import ZptGenWeight
from ZZGenWeight import ZZGenWeight

from DevTools.Utilities.utilities import getCMSSWVersion
from Candidates import *

try:
    from progressbar import ProgressBar, ETA, Percentage, Bar, SimpleProgress
    hasProgress = False
except:
    hasProgress = False

class TreeUpdaterBase(object):
    '''
    Update a tree
    '''

    def __init__(self,**kwargs):
        inputFileNames = kwargs.pop('inputFileNames',[])
        inputTreeDirectory = kwargs.pop('inputTreeDirectory','')
        inputTreeName = kwargs.pop('inputTreeName','AnalysisTree')
        outputFileName = kwargs.pop('outputFileName','analysisTree.root')
        outputTreeName = kwargs.pop('outputTreeName','AnalysisTree')
        self.outputTreeName = outputTreeName
        if hasProgress:
            self.pbar = kwargs.pop('progressbar',ProgressBar(widgets=['{0}: '.format(outputTreeName),' ',SimpleProgress(),' events ',Percentage(),' ',Bar(),' ',ETA()]))
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
        self.treename = '{0}/{1}'.format(inputTreeDirectory,inputTreeName) if inputTreeDirectory else inputTreeName
        self.totalEntries = 0
        self.numLumis = 0
        self.numEvents = 0
        self.summedWeights = 0
        logging.info('Getting information')
        if len(self.fileNames)==0: logging.warning('No files to process')
        if len(self.fileNames)>1: logging.warning('More than one file requested, only processing the first file')
        for f,fName in enumerate(self.fileNames):
            if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
            tfile = ROOT.TFile.Open(fName)
            tree = tfile.Get(self.treename)
            self.totalEntries += tree.GetEntries()
            if not hasattr(self,'version'):
                tree.GetEntry(1)
                if hasattr(tree,'provenance'):
                    ver = tree.provenance[0].split('_')
                    self.version = ''.join([ver[1],ver[2],'X'])
                else:
                    self.version = getCMSSWVersion()
            tfile.Close('R')
        logging.info('Analysis is running with version {0}'.format(self.version))
        self.flush()
        if not len(self.fileNames): raise Exception
        # other input files
        self.pileupWeights = PileupWeights(self.version)
        self.fakeRates = FakeRates(self.version)
        self.leptonScales = LeptonScales(self.version)
        self.triggerScales = TriggerScales(self.version)
        self.triggerPrescales = TriggerPrescales(self.version)
        self.zptGenWeight = ZptGenWeight(self.version)
        self.zzGenWeight = ZZGenWeight(self.version)
        # tfile
        fName = self.fileNames[0]
        if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs' if self.hasHDFS else 'root://cmsxrootd.hep.wisc.edu/',fName)
        self.tfile = ROOT.TFile.Open(fName,'READ')
        self.oldtree = self.tfile.Get(self.treename)
        self.outfile = ROOT.TFile(outputFileName,"recreate")
        self.tree = self.oldtree.CloneTree(0)
        summedWeights = self.tfile.Get('summedWeights')
        self.summedWeights = summedWeights.GetBinContent(1)

    def __exit__(self, type, value, traceback):
        self.finish()

    def __del__(self):
        self.finish()

    def finish(self):
        print ''
        logging.info('Finishing')
        if hasattr(self,'outfile'):
            self.outfile.cd()
            cutflowHist = ROOT.TH1F('summedWeights','summedWeights',1,0,1)
            cutflowHist.SetBinContent(1,self.summedWeights)
            self.outfile.Write()
            self.outfile.Close()
        if hasattr(self,'leptonScales'): self.leptonScales.finish()

    def flush(self):
        sys.stdout.flush()
        sys.stderr.flush()

    ############################
    ### primary updater loop ###
    ############################
    def update(self):
        logging.info('Beginning Updater')
        start = time.time()
        new = start
        old = start
        if hasProgress:
            self.pbar.maxval = self.totalEntries
            self.pbar.start()
            total = 0
            for row in self.oldtree:
                total += 1
                self.pbar.update(total)
                self.perRowAction()
            self.tfile.Close('R')
            self.pbar.update(self.totalEntries)
        else:
            total = 0
            for row in self.oldtree:
                total += 1
                if total==2: start = time.time() # just ignore first event for timing
                if total % 10000 == 1:
                    cur = time.time()
                    elapsed = cur-start
                    remaining = float(elapsed)/total * float(self.totalEntries) - float(elapsed)
                    mins, secs = divmod(int(remaining),60)
                    hours, mins = divmod(mins,60)
                    logging.info('{0}: Processing event {1}/{2} - {3}:{4:02d}:{5:02d} remaining'.format(self.outputTreeName,total,self.totalEntries,hours,mins,secs))
                    self.flush()
                self.perRowAction()
            self.tfile.Close('R')

    def perRowAction(self):
        '''Per row action, can be overridden'''
        pass
        #self.tree.fill()







