#!/usr/bin/env python
import sys
import os

import ROOT

from DevTools.Utilities.utilities import getCMSSWVersion

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    histName = 'pileup'
    fileName = 'pileup/pileup.root'
    
    # read data
    dataFileName = 'pileup/PileUpData.root'
    datafile = ROOT.TFile(dataFileName)
    histdata = datafile.Get(histName)
    
    # now get values
    numbins = histdata.GetNbinsX()
    vals = [0]*numbins
    for b in range(numbins):
        d = histdata.GetBinContent(b+1)
        vals[b] = d

    print vals
    


if __name__ == "__main__":
    status = main()
    sys.exit(status)
                    
