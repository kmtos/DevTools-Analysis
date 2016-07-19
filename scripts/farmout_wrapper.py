#!/usr/bin/env python
import argparse
import sys
import subprocess
import os

# import run script
from DevTools.Analyzer.WZAnalysis import main as runWZ
from DevTools.Analyzer.ZZAnalysis import main as runZZ
from DevTools.Analyzer.DYAnalysis import main as runDY
from DevTools.Analyzer.ZFakeRateAnalysis import main as runZFakeRate
from DevTools.Analyzer.ChargeAnalysis import main as runCharge
from DevTools.Analyzer.TauChargeAnalysis import main as runTauCharge
from DevTools.Analyzer.Hpp3lAnalysis import main as runHpp3l
from DevTools.Analyzer.Hpp4lAnalysis import main as runHpp4l
from DevTools.Analyzer.DijetFakeRateAnalysis import main as runDijetFakeRate
from DevTools.Analyzer.WTauFakeRateAnalysis import main as runWTauFakeRate
from DevTools.Analyzer.WFakeRateAnalysis import main as runWFakeRate
from DevTools.Analyzer.ElectronAnalysis import main as runElectron
from DevTools.Analyzer.MuonAnalysis import main as runMuon
from DevTools.Analyzer.TauAnalysis import main as runTau

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit analyzers')

    parser.add_argument('analysis', type=str, choices=['ZZ', 'WZ', 'DY', 'ZFakeRate', 'Charge', 'TauCharge', 'Hpp3l', 'Hpp4l', 'Electron', 'Muon','Tau', 'DijetFakeRate', 'WTauFakeRate', 'WFakeRate'], help='Analysis to submit')

    return parser.parse_args(argv)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # add on the input and output
    argv = ['--inputFileList',os.environ['INPUT'],'--outputFile',os.environ['OUTPUT']]

    # run the analyzer
    if args.analysis=='WZ':
        status = runWZ(argv)
    elif args.analysis=='ZZ':
        status = runZZ(argv)
    elif args.analysis=='DY':
        status = runDY(argv)
    elif args.analysis=='ZFakeRate':
        status = runZFakeRate(argv)
    elif args.analysis=='Charge':
        status = runCharge(argv)
    elif args.analysis=='TauCharge':
        status = runTauCharge(argv)
    elif args.analysis=='Hpp3l':
        status = runHpp3l(argv)
    elif args.analysis=='Hpp4l':
        status = runHpp4l(argv)
    elif args.analysis=='DijetFakeRate':
        status = runDijetFakeRate(argv)
    elif args.analysis=='WTauFakeRate':
        status = runWTauFakeRate(argv)
    elif args.analysis=='WFakeRate':
        status = runWFakeRate(argv)
    elif args.analysis=='Electron':
        status = runElectron(argv)
    elif args.analysis=='Muon':
        status = runMuon(argv)
    elif args.analysis=='Tau':
        status = runTau(argv)
    else:
        status = 0

    return status

if __name__ == "__main__":
    status = main()
    sys.exit(status)
