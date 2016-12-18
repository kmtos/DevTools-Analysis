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
from DevTools.Analyzer.TriggerCountAnalysis import main as runTriggerCount

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit analyzers')

    parser.add_argument('analysis', type=str, choices=['ZZ', 'WZ', 'DY', 'ZFakeRate', 'Charge', 'TauCharge', 'Hpp3l', 'Hpp4l', 'Electron', 'Muon','Tau', 'DijetFakeRate', 'WTauFakeRate', 'WFakeRate','TriggerCount'], help='Analysis to submit')
    parser.add_argument('shift', type=str, nargs='?', default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')

    return parser.parse_args(argv)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # add on the input and output
    argv = ['--inputFileList',os.environ['INPUT'],'--outputFile',os.environ['OUTPUT'],'--shift',args.shift]

    # run the analyzer
    if args.analysis=='WZ':
        func = runWZ
    elif args.analysis=='ZZ':
        func = runZZ
    elif args.analysis=='DY':
        func = runDY
    elif args.analysis=='ZFakeRate':
        func = runZFakeRate
    elif args.analysis=='Charge':
        func = runCharge
    elif args.analysis=='TauCharge':
        func = runTauCharge
    elif args.analysis=='Hpp3l':
        func = runHpp3l
    elif args.analysis=='Hpp4l':
        func = runHpp4l
    elif args.analysis=='DijetFakeRate':
        func = runDijetFakeRate
    elif args.analysis=='WTauFakeRate':
        func = runWTauFakeRate
    elif args.analysis=='WFakeRate':
        func = runWFakeRate
    elif args.analysis=='Electron':
        func = runElectron
    elif args.analysis=='Muon':
        func = runMuon
    elif args.analysis=='Tau':
        func = runTau
    elif args.analysis=='TriggerCount':
        func = runTriggerCount
    else:
        return 0

    status = func(argv)

    return status

if __name__ == "__main__":
    func = main()
    sys.exit(func)
