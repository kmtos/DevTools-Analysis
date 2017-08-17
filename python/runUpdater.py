#!/usr/bin/env python

# import run script
from DevTools.Analyzer.EGTreeUpdater import main as runEG
from DevTools.Analyzer.ThreePhotonTreeUpdater import main as runThreePhoton

funcMap = {
    'EG'         : runEG,
    'ThreePhoton': runThreePhoton,
}

def runUpdater(analysis,argv):
    '''Return analysis function'''
    if analysis in funcMap:
        return funcMap[analysis](argv)
    else:
        return 0
