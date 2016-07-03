# common utilities for analyzers
import os
import sys
import glob

import ROOT

from DevTools.Utilities.utilities import ZMASS, getCMSSWVersion

def deltaPhi(phi0,phi1):
    result = phi0-phi1
    while result>ROOT.TMath.Pi():
        result -= 2*ROOT.TMath.Pi()
    while result<=-ROOT.TMath.Pi():
        result += 2*ROOT.TMath.Pi()
    return result

def deltaR(eta0,phi0,eta1,phi1):
    deta = eta0-eta1
    dphi = deltaPhi(phi0,phi1)
    return ROOT.TMath.Sqrt(deta**2+dphi**2)

latestNtuples = {
    '76X' : '2016-04-22_DevTools_76X_v1',
    '80X' : '2016-06-22_DevTools_80X_v1',
}

def getNtupleDirectory(version=None):
    baseDir = '/hdfs/store/user/dntaylor'
    if not version: version = getCMSSWVersion()
    if version in latestNtuples:
        return os.path.join(baseDir,latestNtuples[version])

def getTestFiles(sample,n=1):

    sampleMap = {
        'wz'  : 'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
        'data': 'DoubleMuon',
        #'hpp' : 'HPlusPlusHMinusMinusHTo4L_M-500_13TeV-pythia8',
        'hpp' : 'HPlusPlusHMinusMinusHTo4L_M-500_TuneCUETP8M1_13TeV_pythia8',
        #'dy'  : 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
        'dy'  : 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'w'   : 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
        'SingleMuon': 'SingleMuon',
    }

    if sample not in sampleMap: return []
    
    files = [f.replace('/hdfs','') for f in glob.glob('{0}/{1}/*/*/*/*.root'.format(getNtupleDirectory(),sampleMap[sample]))]

    return files[:min(n,len(files))]

