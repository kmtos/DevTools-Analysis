import argparse
import logging
import sys

import csv
from Candidates import *
import KinematicFitter

import sys
import itertools
import operator

import ROOT

from DataFormats.FWLite import Handle, Events
from DevTools.Utilities.utilities import *

DeltaR = ROOT.Math.VectorUtil.DeltaR

mu12, mu12Label = Handle("std::vector<pat::Muon>"), "Mu1Mu2"
mu3,  mu3Label = Handle("std::vector<pat::Muon>"), "Mu3"
mets, metsLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
taus, tausLabel = Handle("std::vector<pat::Tau>"), "muHadTauDMIsoSelector"
jets, jetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
genParticles, genParticlesLabel = Handle("std::vector<pat::PackedCandidate>"), "packedPFCandidates"
outputFile  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/BSUB/IRNAME_Plots.root", "recreate")
mumu_Mass   = ROOT.TH1F("mumu_Mass", "", 600, 0, 30)
tautau_Mass = ROOT.TH1F("tautau_Mass", "", 600, 0, 30)
h_Mass      = ROOT.TH1F("h_Mass", "", 3000, 0, 1000)
mumu_Mass_Before    = ROOT.TH1F("mumu_Mass_Before", "", 3000, 0, 30)
DiTauMassVSMVA = ROOT.TH2F("DiTauMassVSMVA", "", 25, 0, 25, 20, -10, 1)
DiTauMassVSTauMu3dR = ROOT.TH2F("DiTauMassVSTauMu3dR", "", 25, 0, 25, 20, 0, 4)
DiTauMassVSDiMudR = ROOT.TH2F("DiTauMassVSDiMudR", "", 25, 0, 25, 20, 0, 4)
DiTauMassVSTauHMu1dR = ROOT.TH2F("DiTauMassVSTauHMu1dR", "", 25, 0, 25, 20, 0, 4)
DiTauMassVSmu3Pt = ROOT.TH2F("DiTauMassVSmu3Pt", "", 25, 0, 25, 20, 0, 500)
DiTauMassVStauPt = ROOT.TH2F("DiTauMassVStauPt", "", 25, 0, 25, 20, 0, 500)
DiTauMassVSDiTauPt = ROOT.TH2F("DiTauMassVSDiTauPt", "", 25, 0, 25, 20, 0, 500)
DiTauMassVSDiMuMass = ROOT.TH2F("DiTauMassVSDiMuMass", "", 25, 0, 25, 20, 0, 25)
NEvents = ROOT.TH1F("NEvents", "",10,-0.5,9.5)
NEvents.GetXaxis().SetBinLabel(1, "All")
NEvents.GetXaxis().SetBinLabel(2, "mu1 pt > 26")
NEvents.GetXaxis().SetBinLabel(1, "Found ditau")
NEvents.GetXaxis().SetBinLabel(1, "All")
NEvents.GetXaxis().SetBinLabel(1, "All")
NEvents.GetXaxis().SetBinLabel(1, "All")
NEvents.GetXaxis().SetBinLabel(1, "All")


 
ROOT.gSystem.Load('libCondFormatsBTauObjects')
ROOT.gSystem.Load('libCondToolsBTau')
calib = ROOT.BTagCalibration("csv", "/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CSVv2_Moriond17_B_H.csv")
reader = ROOT.BTagCalibrationReader(1, "central")  
reader.load(calib, 0, "comb")
reader.load(calib, 1, "comb")
reader.load(calib, 2, "incl")

FILES = ["root://eoscms/eos/cms/store/group/phys_higgs/HiggsExo/ktos/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/NoMuonClean_SIG_h125a5_MuonSelectionOnly_SEP24/180925_021145/0000/egionB_selection_3.root"]
for FILE in FILES:
  for iev,event in enumerate(FILE):
    NEvents.Fill(0)
    print type(mu12), type(mu12Label), type(event)
    event.getByLabel(mu12Label, mu12)
    event.getByLabel(mu3Label, mu3)
    event.getByLabel(metsLabel, mets)
    event.getByLabel(tausLabel, taus)
    event.getByLabel(jetsLabel, jets)
    event.getByLabel(genParticles, genParticlesLabel)

    highestPt = 0
    for i,mu in enumerate(mu12.product()):
      if highestPt == 0: 
        m1 = mu
        highestPt = m1.pt()
      elif mu.pt() > m1.pt():
        m2 = m1
        m1 = mu
        highestPt = m1.pt()
      else:
        m2 = mu
    if m1.pt() < 26.0: 
      continue
    NEvents.Fill(1)
     
    #################################
    ## Check the Gen particle status
    #################################
    bestdR1 = 1000
    bestdR2 = 1000
    foundGenTau = False 
    foundMu3 = False
    for i, genP in enumerate(genParticles.product() ):
      if fabs(genP.pdgId()) == 13:
        dR_1 = DeltaR(m1.p4(), genP.p4() )
        dR_2 = DeltaR(m2.p4(), genP.p4() )
        pt_diff1 = (genP.pt() - m1.pt() ) / genP.pt()
        pt_diff2 = (genP.pt() - m2.pt() ) / genP.pt()
        if dR_1 < bestdR1 and dR_1 < 0.2 and pt_diff1 < 0.1: 
          gen_m1 = genP
          bestdR1 = dR_1

        if dR_2 < bestdR2 and dR_2 < 0.2 and pt_diff2 < 0.1:
          gen_m2 = genP
          bestdR2 = dR_2

        if dR_2 > 0.35 and dR_1 > 0.35 and genP.pt() > 3.0 and  en_P.eta() < 2.4 and gen_P.eta() > -2.4: foundMu3 = True

      if fabs(genP.pdgId()) == 15 and not foundGenTau:
        for i in genP.daughters():
          print "\tDAUGHTER: pdgId=", i.pdgId(), "  pt=", i.pt(), "  eta=", i.eta()
        if genP.pt() > 20.0 and gen_P.eta() < 2.4 and gen_P.eta() > -2.4: foundGenTau = True
  
 
    if gen_m1.pt() < 26.0 or gen_m1.eta() < -2.4 or gen_m1.eta() > 2.4: continue
    if gen_m2.pt() < 3.0 or gen_m2.eta() < -2.4 or gen_m2.eta() > 2.4: continue

    totalWeight = 1.0 
    if not foundGenTau or not foundMu3: continue
    mumu_Mass_Before.Fill(diMu.M(), totalWeight )

    bestdR = 10000
    foundTau = False
    for i, tau in enumerate(taus.product() ):
      if tau.pt() < 20.0 or tau.eta() > 2.4 or tau.eta() < -2.4 or tau.ID("decayModeFinding") < 0.5: continue
      dr1 = DeltaR(tau.p4(), m1.p4() )
      dr2 = DeltaR(tau.p4(), m2.p4() )
      if dr1 < 0.8 or dr2 < 0.8: continue
      for i, mu in enumerate(mu3.product() ):
        dr3 = DeltaR(tau.p4(), mu.p4() )
        dr13 = DeltaR(mu.p4(), m1.p4() )
        dr23 = DeltaR(mu.p4(), m2.p4() )
        if dr3 < bestdR and dr3 < 0.8 and dr13 > 0.4 and dr23 > 0.4:
          foundTau = True
          bestdR = dr3
          tm = mu
          th = tau
          mva = tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT")
          minmva = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw")
      bestTauJetdR = 100
      for i, jet in enumerate(jets.product() ):
        drjet = DeltaR(tau.p4(), jet.p4() )
        if drjet < bestTauJetdR:
          bestTauJetdR = drjet
          tj = jet

    if not foundTau: continue
    if mva > 0.5 and minmva > -0.5 : continue

    flavor = tj.hadronFlavour()
    if flavor == 5:   csvSF = reader.eval(0, tj.eta(), tj.pt() )
    elif flavor == 4: csvSF = reader.eval(1, tj.eta(), tj.pt() )
    else:             csvSF = reader.eval(2, tj.eta(), tj.pt() )
    if tj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > .8484: continue

    NEvents.Fill(2)
    diTau = ROOT.TLorentzVector()
    diMu = ROOT.TLorentzVector()
    diMu = m1.p4() + m2.p4()
    diTau = tm.p4() + th.p4()

    fourBody = m1.p4() + m2.p4() + tm.p4() + th.p4()
    mumu_Mass.Fill(   diMu.M(), totalWeight )
    tautau_Mass.Fill( diTau.M(), totalWeight )
    h_Mass.Fill(      fourBody.M(), totalWeight )
    DiTauMassVSMVA.Fill(diTau.M(), mva)
    DiTauMassVSTauMu3dR.Fill(diTau.M(), bestdR)
    DiTauMassVSDiMudR.Fill(diTau.M(), DeltaR(m1.p4(), m2.p4()) )
    DiTauMassVSTauHMu1dR.Fill(diTau.M(), DeltaR(m1.p4(),th.p4()) )
    DiTauMassVSmu3Pt.Fill(diTau.M(), tm.pt())
    DiTauMassVStauPt.Fill(diTau.M(), th.pt())
    DiTauMassVSDiTauPt.Fill(diTau.M(), diTau.Pt())
    DiTauMassVSDiMuMass.Fill(diTau.M(), diMu.M())


outputFile.cd()
mumu_Mass.Write()
mumu_Mass_Before.Write() 
tautau_Mass.Write()
h_Mass.Write()
TotalWeight.Write()
DiTauMassVSMVA.Write()
DiTauMassVSTauMu3dR.Write()
DiTauMassVSDiMudR.Write()
DiTauMassVSTauHMu1dR.Write()
DiTauMassVSmu3Pt.Write()
DiTauMassVStauPt.Write()
DiTauMassVSDiTauPt.Write()
DiTauMassVSDiMuMass.Write()
outputFile.Write()
outputFile.Close()
csvFile.close()
print "Ratio of tau selection to no tau selection=", mumu_Mass.Integral() / mumu_Mass_Before.Integral()
