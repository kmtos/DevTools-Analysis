#!/usr/bin/env python

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
from DevTools.Analyzer.SingleJetAnalysis import main as runSingleJet


def runAnalysis(analysis,argv):
    '''Return analysis function'''
    if analysis=='WZ':
        func = runWZ
    elif analysis=='ZZ':
        func = runZZ
    elif analysis=='DY':
        func = runDY
    elif analysis=='ZFakeRate':
        func = runZFakeRate
    elif analysis=='Charge':
        func = runCharge
    elif analysis=='TauCharge':
        func = runTauCharge
    elif analysis=='Hpp3l':
        func = runHpp3l
    elif analysis=='Hpp4l':
        func = runHpp4l
    elif analysis=='DijetFakeRate':
        func = runDijetFakeRate
    elif analysis=='WTauFakeRate':
        func = runWTauFakeRate
    elif analysis=='WFakeRate':
        func = runWFakeRate
    elif analysis=='Electron':
        func = runElectron
    elif analysis=='Muon':
        func = runMuon
    elif analysis=='Tau':
        func = runTau
    elif analysis=='TriggerCount':
        func = runTriggerCount
    elif analysis=='ThreeLepton':
        func = runThreeLepton
    elif analysis=='FourPhoton':
        func = runFourPhoton
    elif analysis=='ThreePhoton':
        func = runThreePhoton
    elif analysis=='TwoPhoton':
        func = runTwoPhoton
    elif analysis=='EG':
        func = runEG
    elif analysis=='DYGG':
        func = runDYGG
    elif analysis=='MMG':
        func = runMMG
    elif analysis=='SingleJet':
        func = runSingleJet
    else:
        return 0

    return func(argv)

