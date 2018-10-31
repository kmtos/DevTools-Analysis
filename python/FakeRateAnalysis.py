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

def isAncestor(a,p) :
  if a == p : 
    return True
  for i in xrange(0,p.numberOfMothers()) :
    if isAncestor(a,p.mother(i)) :
      return True
  return False


DeltaR = ROOT.Math.VectorUtil.DeltaR
ENDFILENUM = ENDFILENUMVALUE

mu12, mu12Label = Handle("std::vector<pat::Muon>"), "Mu1Mu2"
mu3,  mu3Label = Handle("std::vector<pat::Muon>"), "Mu3"
mets, metsLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
if "NoMuonClean" in "FILE_PATH": 
  print "noMuonClean"
  taus, tausLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
else:  taus, tausLabel = Handle("std::vector<pat::Tau>"), "slimmedTausMuonCleaned" 
jets, jetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
prunedParticles, prunedParticlesLabel = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"
packedParticles, packedParticlesLabel = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
outputFile  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/BSUB/DIRNAME/DIRNAME_Plots.root", "recreate")
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
FailBefore = ROOT.TH1F("FailBefore", "",10,-0.5,9.5)
FailBefore.GetXaxis().SetBinLabel(1, "All")
FailBefore.GetXaxis().SetBinLabel(2, "mu1")
FailBefore.GetXaxis().SetBinLabel(3, "mu2")
FailBefore.GetXaxis().SetBinLabel(4, "mu3")
FailBefore.GetXaxis().SetBinLabel(5, "tau")
FailBefore.GetXaxis().SetBinLabel(6, "nMuons")
FailAfter = ROOT.TH1F("FailAfter", "",10,-0.5,9.5)
FailAfter.GetXaxis().SetBinLabel(1, "All Events")
FailAfter.GetXaxis().SetBinLabel(2, "AllTaus")
FailAfter.GetXaxis().SetBinLabel(3, "pt/Eta")
FailAfter.GetXaxis().SetBinLabel(4, "DM")
FailAfter.GetXaxis().SetBinLabel(5, "NoMu3")
FailAfter.GetXaxis().SetBinLabel(6, "TooClose")
FailAfter.GetXaxis().SetBinLabel(7, "No Mu3")
FailAfter.GetXaxis().SetBinLabel(8, "MVA")
FailAfter.GetXaxis().SetBinLabel(9, "PassMVA")


 
ROOT.gSystem.Load('libCondFormatsBTauObjects')
ROOT.gSystem.Load('libCondToolsBTau')
calib = ROOT.BTagCalibration("csv", "/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CSVv2_Moriond17_B_H.csv")
reader = ROOT.BTagCalibrationReader(1, "central")  
reader.load(calib, 0, "comb")
reader.load(calib, 1, "comb")
reader.load(calib, 2, "incl")

FILES = []
 
for i in range(1,ENDFILENUM):
  if ROOT.gSystem.AccessPathName("FILE_PATHRegionB_selection_"+str(i)+".root") == 0:
    FILES.append( Events("root://eoscms/FILE_PATHRegionB_selection_"+str(i)+".root"))

print "NUM FILES=", len(FILES)
for FILE in FILES:
  for iev,event in enumerate(FILE):
    FailBefore.Fill(0)
    event.getByLabel(mu12Label, mu12)
    event.getByLabel(mu3Label, mu3)
    event.getByLabel(metsLabel, mets)
    event.getByLabel(tausLabel, taus)
    event.getByLabel(jetsLabel, jets)
    event.getByLabel(prunedParticlesLabel, prunedParticles)
    event.getByLabel(packedParticlesLabel, packedParticles)

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
    diMu = ROOT.TLorentzVector()
    diMu = m1.p4() + m2.p4()
     
    ###################################
    # Get Pruned Higgs  Particle
    ###################################
    print "\n\nNewEvent"
    for i, genP in enumerate(prunedParticles.product() ):
      if (genP.pdgId() == 25 or genP.pdgId() == 35):
        higgs = genP
        break

    ###################################
    # Get Finanl Pruned Higgs Particle
    ###################################
    foundDaughter = True
    while foundDaughter:
      foundDaughter = False
      for i in range(higgs.numberOfDaughters() ):
        print "Higgs Duaghter #" + str(i) + "=", higgs.daughter(i).pdgId(), higgs.daughter(i).pt(), higgs.daughter(i).eta()
        if abs(higgs.daughter(i).pdgId() ) == 25 or abs(higgs.daughter(i).pdgId() ) == 35:
          foundDaughter = True
          higgs = higgs.daughter(i) 
          break

    a1s = []
    for i in range(higgs.numberOfDaughters() ):
      if abs(higgs.daughter(i).pdgId() ) == 36: a1s.append(higgs.daughter(i))
    if len(a1s) != 2:
      for i in range(higgs.numberOfDaughters() ): print "FAIL: Higgs daughter #" + str(i) + "=", higgs.daughter(i).pdgId(), higgs.daughter(i).pt(), higgs.daughter(i).eta()
      break 

    ###################################
    # Get Final Pruned a1 Particle
    ###################################
    for a1_ite in range(len(a1s) ):
      foundDaughter = True
      while foundDaughter:
        foundDaughter = False
        for i in range(a1s[a1_ite].numberOfDaughters() ):
          print "a1 #" + str(a1_ite) +" Daughter #" + str(i) + "=", a1s[a1_ite].daughter(i).pdgId(), a1s[a1_ite].daughter(i).pt(), a1s[a1_ite].daughter(i).eta()
          if abs(a1s[a1_ite].daughter(i).pdgId() ) == 36:
            foundDaughter = True
            a1s[a1_ite] = a1s[a1_ite].daughter(i)
            break
          
    a1muons = []
    a1taus  = []
    for a1 in a1s:
      for i in range(a1.numberOfDaughters() ):
        if abs(a1.daughter(i).pdgId() ) == 13: 
          a1muons.append(a1.daughter(i))
        if abs(a1.daughter(i).pdgId() ) == 15:
          a1taus.append(a1.daughter(i))
    if len(a1muons) != 2 or len(a1taus) != 2:
      for i in range(a1[0].numberOfDaughters() ): print "FAIL: a1[0] Duaghter #" + str(i) + "=", a1[0].daughter(i).pdgId(), a1[0].daughter(i).pt(), a1[0].daughter(i).eta()
      for i in range(a1[1].numberOfDaughters() ): print "FAIL: a1[1] Duaghter #" + str(i) + "=", a1[1].daughter(i).pdgId(), a1[1].daughter(i).pt(), a1[1].daughter(i).eta()
      break      
 
    #######################
    # Get Tau Constituents
    #######################
    for a1tau_ite in range(len(a1taus) ):
      foundDaughter = True
      while foundDaughter:
        foundDaughter = False
        for i in range(a1taus[a1tau_ite].numberOfDaughters() ):
          print "Tau #" + str(a1tau_ite) + " Daughter #" + str(i) + "=", a1taus[a1tau_ite].daughter(i).pdgId(), a1taus[a1tau_ite].daughter(i).pt(), a1taus[a1tau_ite].daughter(i).eta()
          if abs(a1taus[a1tau_ite].daughter(i).pdgId() ) == 15:
            foundDaughter = True
            a1taus[a1tau_ite] = a1taus[a1tau_ite].daughter(i)
            break

    foundTauMuPtEta = False
    foundTauHadPtEta = False
    foundTauElectron = False
    for tau in range(len(a1taus)):
      tauConst = ROOT.TLorentzVector()
      foundMuonDecay = False
      addedP4 = False
      for i in range(a1taus[tau].numberOfDaughters() ):
        daughter = a1taus[tau].daughter(i)
        print "Constituents: Tau #" + str(tau) + " daughter #" + str(i) + " pdgId()=", daughter.pdgId()

        if abs(daughter.pdgId() ) == 12 or abs(daughter.pdgId() ) == 14 or abs(daughter.pdgId() ) == 16: continue

        elif abs(daughter.pdgId() ) == 11:
          print "FAIL Electronic  Decay"
          foundTauElectron = True
          break

        elif abs(daughter.pdgId() ) == 13:
          foundMuonDecay = True
          tau_mu = daughter
          break

        elif not addedP4:
          tauConst  = daughter.p4()
          addedP4 = True

        else:
          tauConst += daughter.p4()

      if foundMuonDecay: 
        if tau_mu.pt() > 3.0 and abs(tau_mu.eta() ) < 2.4:
          foundTauMuPtEta = True
      else:
        print "Total pt=", tauConst.Pt(), " eta=", tauConst.Eta()
        if tauConst.Pt() > 20.0 and abs(tauConst.Eta()) < 2.4: foundTauHadPtEta = True
#        if "NoMuonClean" in "FILE_PATH" and tauConst.Pt() > 20.0 and abs(tauConst.Eta()) < 2.4: foundTauHadPtEta = True
#        elif "NoMuonClean" not in "FILE_PATH" and tauConst.Pt() > 10.0 and abs(tauConst.Eta()) < 2.4: foundTauHadPtEta = True

    print "RESULTS"
    if foundTauElectron: continue

    if a1muons[0].pt() < a1muons[1].pt(): 
      temp = a1muons[0]
      a1muons[0] = a1muons[1]
      a1muons[1] = temp
    for i in range(len(a1muons)): print "MUON #" + str(i) + " pt=", a1muons[i].pt(),  "\teta=", a1muons[i].eta()
    
    if a1muons[0].pt() <  26.0 and abs(a1muons[0].eta()) > 2.4: 
      print "FAIL Mu1"
      FailBefore.Fill(1)
      continue
    elif a1muons[1].pt() <  3.0 and abs(a1muons[1].eta()) > 2.4: 
      print "FAIL Mu2"
      FailBefore.Fill(2)
      continue
    elif not foundTauMuPtEta:
      print "FAIL Mu3"
      FailBefore.Fill(3)
      continue
    elif not foundTauHadPtEta: 
      print "FAIL tau"
      FailBefore.Fill(4)
      continue
    else:
      mumu_Mass_Before.Fill(diMu.M())
       
  

    #########################
    #Time to check numerator
    #########################
    bestdR = 10000
    foundTau = False
    FailAfter.Fill(0)
    for i, tau in enumerate(taus.product() ):
      FailAfter.Fill(1)

#      if "NoMuonClean" in "FILE_PATH":
      if tau.pt() < 20.0 or tau.eta() > 2.4 or tau.eta() < -2.4:
        FailAfter.Fill(2)
        continue

#      if "NoMuonClean" not in "FILE_PATH":
#        if tau.pt() < 10.0 or tau.eta() > 2.4 or tau.eta() < -2.4:
#          FailAfter.Fill(2)
#          continue


      if tau.tauID("decayModeFinding") < 0.5: 
        FailAfter.Fill(3)
        continue

      dr1 = DeltaR(tau.p4(), m1.p4() )
      dr2 = DeltaR(tau.p4(), m2.p4() )
      if dr1 < 0.8 or dr2 < 0.8: 
        FailAfter.Fill(4)
        continue

      for i, mu in enumerate(mu3.product() ):
        dr3 = DeltaR(tau.p4(), mu.p4() )
        if dr3 < bestdR and dr3 < 0.8:
          foundTau = True
          bestdR = dr3
          tm = mu
          th = tau
          mva = tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT")
          minmva = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw")
      bestTauJetdR = 100
#      for i, jet in enumerate(jets.product() ):
#        drjet = DeltaR(tau.p4(), jet.p4() )
#        if drjet < bestTauJetdR:
#          bestTauJetdR = drjet
#          tj = jet

    if not foundTau: 
      FailAfter.Fill(5)
      continue
    if mva < 0.5 and minmva < -0.5 : 
      FailAfter.Fill(6)
      continue

#    flavor = tj.hadronFlavour()
#    if flavor == 5:   csvSF = reader.eval(0, tj.eta(), tj.pt() )
#    elif flavor == 4: csvSF = reader.eval(1, tj.eta(), tj.pt() )
#    else:             csvSF = reader.eval(2, tj.eta(), tj.pt() )
#    if tj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > .8484: 
#      FailAfter.Fill(7)
#      continue
#    FailAfter.Fill(8)

    diTau = ROOT.TLorentzVector()
    diTau = tm.p4() + th.p4()

    fourBody = m1.p4() + m2.p4() + tm.p4() + th.p4()
    mumu_Mass.Fill(   diMu.M())
    tautau_Mass.Fill( diTau.M())
    h_Mass.Fill(      fourBody.M())
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
print "Ratio of tau selection to no tau selection=", mumu_Mass.Integral() / mumu_Mass_Before.Integral(), "\t", mumu_Mass_Before.Integral()
tautau_Mass.Write()
h_Mass.Write()
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
