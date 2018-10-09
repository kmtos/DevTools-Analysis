import os

from DevTools.Utilities.hdfsUtils import *
from DevTools.Utilities.utilities import *

import ROOT

job = 'crab_2018-10-05_LumiTree_MuMuTauTau_80X_v1'

amasses = ['3p6', 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21]
hmasses = [125, 300, 750]
only125 = ['3p6', 4, 6, 8, 10, 12, 14]

samples = []
for h in hmasses:
    for a in amasses:
        if h==125:
           samples += ['SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'.format(a=a)]
        else:
           if a in only125: continue
           samples += ['SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-{h}_M-{a}_TuneCUETP8M1_13TeV_madgraph_pythia8'.format(h=h,a=a)]

for sample in samples:
    baseDir = '/hdfs/store/user/dntaylor/{}/{}'.format(job,sample)
    files = get_hdfs_root_files(baseDir)

    tree = ROOT.TChain('lumiTree/LumiTree')
    for fName in files:
        if fName.startswith('/store'): fName = '{0}/{1}'.format('/hdfs',fName)
        tree.Add(fName)


    summedWeights = {-1: 0}
    for row in tree:
        nw = len(row.summedGenWeights)
        summedWeights[-1] += row.summedWeights
        for i in range(nw):
            if i not in summedWeights: summedWeights[i] = 0
            summedWeights[i] += row.summedGenWeights[i]

    dumpResults(summedWeights,'MuMuTauTau','{}/summedWeights'.format(sample))
