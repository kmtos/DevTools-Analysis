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
from DevTools.Analyzer.ThreeLeptonAnalysis import main as runThreeLepton
from DevTools.Analyzer.FourPhotonAnalysis import main as runFourPhoton
from DevTools.Analyzer.ThreePhotonAnalysis import main as runThreePhoton
from DevTools.Analyzer.TwoPhotonAnalysis import main as runTwoPhoton
from DevTools.Analyzer.EGAnalysis import main as runEG
from DevTools.Analyzer.DYGGAnalysis import main as runDYGG
from DevTools.Analyzer.MMGAnalysis import main as runMMG

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit analyzers')

    # hack to deal with crab format
    argv = argv[1:] if argv[0].isdigit() else argv # first one is a job num
    argv = ['--'+a if not a.startswith('--') else a for a in argv]

    parser.add_argument('--analysis', type=str, choices=['ZZ', 'WZ', 'DY', 'ZFakeRate', 'Charge', 'TauCharge', 'Hpp3l', 'Hpp4l', 'Electron', 'Muon','Tau', 'DijetFakeRate', 'WTauFakeRate', 'WFakeRate','TriggerCount','ThreeLepton','ThreePhoton','TwoPhoton','EG','MMG','FourPhoton','DYGG'], help='Analysis to submit')
    parser.add_argument('--shift', type=str, nargs='?', default='', choices=['','ElectronEnUp','ElectronEnDown','MuonEnUp','MuonEnDown','TauEnUp','TauEnDown','JetEnUp','JetEnDown','JetResUp','JetResDown','UnclusteredEnUp','UnclusteredEnDown'], help='Energy shift')
    parser.add_argument('--outputFile', type=str, nargs='?', default='crab.root', help='Output filename')

    return parser.parse_args(argv)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # add on the input and output
    # import the crab generated PSet.py
    argv = []
    if os.path.isfile('PSet.py'):
        import PSet
        fileList = list(PSet.process.source.fileNames)

        # add on the inputfiles
        argv += ['--inputFiles'] + fileList
    else:
        print 'Warning: Failed to find PSet.py'

    argv += ['--outputFile',args.outputFile,'--shift',args.shift]

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
    elif args.analysis=='ThreeLepton':
        func = runThreeLepton
    elif args.analysis=='FourPhoton':
        func = runFourPhoton
    elif args.analysis=='ThreePhoton':
        func = runThreePhoton
    elif args.analysis=='TwoPhoton':
        func = runTwoPhoton
    elif args.analysis=='EG':
        func = runEG
    elif args.analysis=='DYGG':
        func = runDYGG
    elif args.analysis=='MMG':
        func = runMMG
    else:
        return 0

    status = func(argv)

    # write xml file
    with open('FrameworkJobReport.xml','w') as f:
        f.write(
'''<FrameworkJobReport>
<ReadBranches>
</ReadBranches>
<PerformanceReport>
  <PerformanceSummary Metric="StorageStatistics">
    <Metric Name="Parameter-untracked-bool-enabled" Value="true"/>
    <Metric Name="Parameter-untracked-bool-stats" Value="true"/>
    <Metric Name="Parameter-untracked-string-cacheHint" Value="application-only"/>
    <Metric Name="Parameter-untracked-string-readHint" Value="auto-detect"/>
    <Metric Name="ROOT-tfile-read-totalMegabytes" Value="0"/>
    <Metric Name="ROOT-tfile-write-totalMegabytes" Value="0"/>
  </PerformanceSummary>
</PerformanceReport>
<GeneratorInfo>
</GeneratorInfo>
</FrameworkJobReport>'''
        )

    return status

if __name__ == "__main__":
    func = main()
    sys.exit(func)
