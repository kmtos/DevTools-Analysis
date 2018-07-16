# common utilities for analyzers
import os
import sys
import glob

import ROOT

from DevTools.Utilities.utilities import ZMASS, getCMSSWVersion
from DevTools.Utilities.hdfsUtils import get_hdfs_root_files

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
    '76X'                  : '2016-09-28_DevTools_76X_v1',
    #'80X'                  : '2016-07-20_DevTools_80X_v1', # ICHEP 2016
    '80X'                  : '2017-11-24_DevTools_80X_v1', # Moriond 2017
    '80XPhoton'            : '2017-11-24_DevTools_80X_v1', # Moriond 2017
    '80XHpp'               : '2018-06-19_DevTools_80X_Hpp_v1',
    #'80XMuMuTauTau'        : '2018-01-31_DevTools_MuMuTauTauSkim_80X_v1', # Moriond 2017
    '80XMuMuTauTau'        : '2018-05-01_DevTools_MuMuTauTauSkim_80X_v1', # Moriond 2017 (add hadronFlavour)
    '80XMuMuTauTauPacked'  : '2018-06-18_DevTools_MuMuTauTauSkim_80X_v2', # Add packed
    '80XMuMu'              : '2018-05-09_DevTools_80X_v1', # Moriond 2017
    #'80XMuMuTauTauUnclean' : '2018-03-09_DevTools_MuMuTauTauSkim_80X_uncleaned_v1', # Moriond 2017
    '80XMuMuTauTauZSkim'   : '2018-01-31_DevTools_MuMuTauTauZSkim_80X_v1', # Moriond 2017
    '94X'                : '2018-03-15_DevTools_94X_v3', # Moriond 2018
}

def getNtupleDirectory(version=None):
    baseDir = '/hdfs/store/user/dntaylor'
    if not version: version = getCMSSWVersion()
    if version in latestNtuples:
        return os.path.join(baseDir,latestNtuples[version])

def getTestFiles(sample,n=1,version=None):
    if not version: version = getCMSSWVersion()

    sampleMap = {
        'wz'    : 'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
        'zz'    : 'ZZTo4L_13TeV_powheg_pythia8',
        'data'  : 'DoubleMuon',
        'hpp'   : 'HPlusPlusHMinusMinusHTo4L_M-500_13TeV-pythia8' if version=='76X' else 'HPlusPlusHMinusMinusHTo4L_M-500_TuneCUETP8M1_13TeV_pythia8',
        'hpp4l' : 'HPlusPlusHMinusMinusHTo4L_M-500_13TeV-pythia8' if version=='76X' else 'HPlusPlusHMinusMinusHTo4L_M-500_TuneCUETP8M1_13TeV_pythia8',
        'hppr4l': 'HPlusPlusHMinusMinusHRTo4L_M-500_13TeV-pythia8' if version=='76X' else 'HPlusPlusHMinusMinusHRTo4L_M-500_TuneCUETP8M1_13TeV-pythia8',
        'hpp3l' : 'HPlusPlusHMinusHTo3L_M-500_TuneCUETP8M1_13TeV_calchep-pythia8' if version=='76X' else 'HPlusPlusHMinusHTo3L_M-500_13TeV-calchep-pythia8',
        'dy'    : 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
        #'dy'    : 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
        'w'     : 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
        'qcd'   : 'QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8',
        'gjet'  : 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8',
        'SingleMuon'    : 'SingleMuon',
        'SingleElectron': 'SingleElectron',
        'DoubleMuon'    : 'DoubleMuon',
        'DoubleEG'      : 'DoubleEG',
        'MuonEG'        : 'MuonEG',
        'Tau'           : 'Tau',
        'haa'           : 'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-15_TuneCUETP8M1_13TeV_madgraph_pythia8',
        'hzz'           : 'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8',
        'hgg'           : 'GluGluHToGG_M-125_13TeV_powheg_pythia8',
    }

    for h in [125,300,750]:
        for a in ['3p6',4,5,6,7,9,11,13,15,17,19,21]:
            if h == 125:
                name = 'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'
            else:
                name = 'SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{h}_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'
            sampleMap['haa_{h}_{a}'.format(h=h,a=a)] = name.format(h=h,a=a)

    if sample not in sampleMap: return []
    
    files = get_hdfs_root_files('{0}/{1}'.format(getNtupleDirectory(version=version),sampleMap[sample]))

    return files[:min(n,len(files))]

