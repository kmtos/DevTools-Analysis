# common utilities for analyzers
import os
import sys
import ROOT

from DevTools.Utilities.utilities import ZMASS

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
    'DevTools' : '2016-04-22_DevTools_76X_v1',
}

def getNtupleDirectory(analysis):
    baseDir = '/hdfs/store/user/dntaylor'
    if analysis in latestNtuples and latestNtuples[analysis]:
        return os.path.join(baseDir,latestNtuples[analysis])

def getTestFiles(type):
    if type=='MC':
        return ['/store/user/dntaylor/2016-04-22_DevTools_76X_v1/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/2016-04-22_DevTools_76X_v1/160422_113101/0000/miniTree_1.root']
    elif type=='Data':
        return ['/store/user/dntaylor/2016-04-22_DevTools_76X_v1/MuonEG/2016-04-22_DevTools_76X_v1/160422_113532/0000/miniTree_1.root']
    elif type=='hpp':
        return ['/store/user/dntaylor/2016-04-22_DevTools_76X_v1/HPlusPlusHMinusMinusHTo4L_M-500_13TeV-pythia8/2016-04-22_DevTools_76X_v1/160422_113847/0000/miniTree_1.root']
    elif type=='dy':
        return ['/store/user/dntaylor/2016-04-22_DevTools_76X_v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016-04-22_DevTools_76X_v1/160422_112025/0000/miniTree_1.root']
    elif type=='long':
        return [
        ]
    else:
        return ['']
