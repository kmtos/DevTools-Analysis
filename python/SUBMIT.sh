#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cd /afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/
cp /afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/BSUB/DIRNAME/ANALYZER.py .
python ANALYZER.py
rm  /afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/ANALYZER.py
exit 0
