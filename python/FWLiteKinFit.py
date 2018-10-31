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

ISMC = True
LUMI = 35900
XSEC = XSECVALUE
SUMMEDWEIGHTS = SUMMEDWEIGHTSVALUE
ENDFILENUM = ENDFILENUMVALUE
APPLYFAKERATE = False

mu12, mu12Label = Handle("std::vector<pat::Muon>"), "Mu1Mu2"
mu3,  mu3Label = Handle("std::vector<pat::Muon>"), "Mu3"
mets, metsLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
taus, tausLabel = Handle("std::vector<pat::Tau>"), "muHadTauDMIsoSelector"
jets, jetsLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
outputFile  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/BSUB/DIRNAME/DIRNAME_Plots.root", "recreate")
mumu_Mass   = ROOT.TH1F("mumu_Mass", "", 600, 0, 30)
tautau_Mass = ROOT.TH1F("tautau_Mass", "", 600, 0, 30)
h_Mass      = ROOT.TH1F("h_Mass", "", 3000, 0, 1000)
hKF_Mass    = ROOT.TH1F("hKF_Mass", "", 3000, 0, 1000)
hKF_ChiSq   = ROOT.TH1F("hKF_ChiSq", "", 30000, 0, 300)
TotalWeight = ROOT.TH1F("TotalWeight", "", 1000, 0, 20)
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


x   = ROOT.RooRealVar("x","x", 0, 30)
y   = ROOT.RooRealVar("y","y", 0, 30)
y1  = ROOT.RooRealVar("y1","y1", 0, 1200)
y2  = ROOT.RooRealVar("y2","y2", 0, 1200)
chi = ROOT.RooRealVar("chi", "chi", 0, 1000)
w   = ROOT.RooRealVar("w","w", 0, 30)
mumumass_dataset         = ROOT.RooDataSet("mumumass_dataset",         "mumumass_dataset", ROOT.RooArgSet(x,y,w) )
mumutautaumass_dataset   = ROOT.RooDataSet("mumutautaumass_dataset",   "mumutautaumass_dataset", ROOT.RooArgSet(x,y,w) )
mumufourBodymass_dataset = ROOT.RooDataSet("mumufourBodymass_dataset", "mumufourBodymass_dataset", ROOT.RooArgSet(x,y1,w) )
mumufourBodyKinFitmass_dataset = ROOT.RooDataSet("mumufourBodyKinFitmass_dataset", "mumufourBodyKinFitmass_dataset", ROOT.RooArgSet(x,y2,chi,w) )
 
if ISMC:
  pileup, pileupLabel = Handle("std::vector<PileupSummaryInfo>"),  "slimmedAddPileupInfo"
  genEvent, genEventLabel = Handle("GenEventInfoProduct"), "generator"
#  PileupFile = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/PileupWeights_69200.root")
#  Pileup_ = PileupFile.Get("PileupWeights")
  PileupFile = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/pileup_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6.root")
  Pileup_ = PileupFile.Get("pileup_scale")

  _fileIDs_BToF  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunBCDEF_ID.root")
  _fileIDs_GH  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunGH_ID.root")
  _fileISOs_BToF  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunBCDEF_ISO.root")
  _fileISOs_GH  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Zmumu_RunGH_ISO.root")
  _fileTrack  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/Tracking_EfficienciesAndSF_BCDEFGH.root")
  _fileTrigger_BToF  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/EfficienciesAndSF_RunBtoF.root")
  _fileTrigger_GH  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/EfficienciesAndSF_Period4.root")
  _file_LowPt  = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/TNP_LowMuPt_EFFICIENCIES.root")
  
  IDsWeight_BToF     =  _fileIDs_BToF.Get("EtavsPtLooseID")
  IDsWeight_GH       =  _fileIDs_GH.Get("EtavsPtLooseID")
  ISOsWeight_BToF    = _fileISOs_BToF.Get("EtavsPtLooseISO")
  ISOsWeight_GH      = _fileISOs_GH.Get("EtavsPtLooseISO")
  TrackWeight_eta    = _fileTrack.Get("ratio_eff_eta3_dr030e030_corr")
  TrackWeight_vtx    = _fileTrack.Get("ratio_eff_vtx_dr030e030_corr")
  TriggerWeight_BToF = _fileTrigger_BToF.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")
  TriggerWeight_GH   = _fileTrigger_GH.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")
  IDsWeight_LowPt_ID = _file_LowPt.Get("hist_EtavsPtLooseID_DatatoMC")
  ISOsWeight_LowPt_ISO = _file_LowPt.Get("hist_EtavsPtLooseISO_DatatoMC")
  ROOT.gSystem.Load('libCondFormatsBTauObjects')
  ROOT.gSystem.Load('libCondToolsBTau')
  calib = ROOT.BTagCalibration("csv", "/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CSVv2_Moriond17_B_H.csv")
  reader = ROOT.BTagCalibrationReader(1, "central")  
  reader.load(calib, 0, "comb")
  reader.load(calib, 1, "comb")
  reader.load(calib, 2, "incl")

#  TrackCorr_eta = []
#  x,y = ROOT.Double(0), ROOT.Double(0)
#  for i in range(TrackWeight_eta.GetN()):
#      TrackWeight_eta.GetPoint(i,x,y)
#      val = {
#          'x'        : float(x),
#          'y'        : float(y),
#          'errx_up'  : float(TrackWeight_eta.GetErrorXhigh(i)),
#          'errx_down': float(TrackWeight_eta.GetErrorXlow(i)),
#          'erry_up'  : float(TrackWeight_eta.GetErrorYhigh(i)),
#          'erry_down': float(TrackWeight_eta.GetErrorYlow(i)),
#      }
#      TrackCorr_eta += [val]
#  for i in range(TrackWeight_vtx.GetN()):
#      TrackWeight_vtx.GetPoint(i,x,y)
#      val = {
#          'x'        : float(x),
#          'y'        : float(y),
#          'errx_up'  : float(TrackWeight_vtx.GetErrorXhigh(i)),
#          'errx_down': float(TrackWeight_vtx.GetErrorXlow(i)),
#          'erry_up'  : float(TrackWeight_vtx.GetErrorYhigh(i)),
#          'erry_down': float(TrackWeight_vtx.GetErrorYlow(i)),
#      }
#      TrackCorr_vtx += [val]
  xTemp = ROOT.Double(0)
  yTemp = ROOT.Double(0)
  TrackCorr_eta = []
  for i in range(TrackWeight_eta.GetN()):
    TrackWeight_eta.GetPoint(i, xTemp, yTemp)
    TrackCorr_eta.append( (float(xTemp), float(yTemp), float(TrackWeight_eta.GetErrorXhigh(i)), float(TrackWeight_eta.GetErrorXlow(i)) ) ) 
  TrackCorr_vtx = []
  for i in range(TrackWeight_vtx.GetN()):
    TrackWeight_vtx.GetPoint(i, xTemp, yTemp)
    TrackCorr_vtx.append( (float(xTemp), float(yTemp), float(TrackWeight_vtx.GetErrorXhigh(i)), float(TrackWeight_vtx.GetErrorXlow(i)) ) )
  print "TrackCorr_eta"
  for i in TrackCorr_eta:
    print i
  print "TrackCorr_vtx"
  for i in TrackCorr_vtx:
    print i
  

elif APPLYFAKERATE:
  fileNameFR = ROOT.TFile("/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/TAUFAKERATES.root")
  histFR_ = fileNameFR.Get("FinalFakeRateDMtoMedIsoOnlyEtavsPt")

FILES = []
 
for i in range(1,ENDFILENUM):
  if ROOT.gSystem.AccessPathName("FILE_PATHRegionB_selection_"+str(i)+".root") == 0:
    FILES.append( Events("root://eoscms/FILE_PATHRegionB_selection_"+str(i)+".root"))

print "NUM FILES=", len(FILES)
csvFile = open('/afs/cern.ch/work/k/ktos/public/DevAnalysis/CMSSW_8_1_0/src/DevTools/Analyzer/python/BSUB/DIRNAME/DIRNAME_EventList.csv','w')
csvFile.write("run,lumi,event,m1Pt,m1Eta,m1Phi,m1Energy,m1Iso,m2Pt,m2Eta,m2Phi,m2Energy,m2Iso,m3Pt,m3Eta,m3Phi,m3Energy,m3Iso,tPt,tEta,tPhi,tEnergy,tMVA,tPassID,tDM,ammMass,ammDeltaR,attMass,attDeltaR,hMass,hMassKinFit,kinFitChi2,genWeight,pileupWeight,triggerWeight,m1IDWeight,m1IsoWeight,m2IDWeight,m2IsoWeight,m3IDWeight,m1TrackWeight,m2TrackWeight,m3TrackWeight,lumiWeight,bVetoWeight,weight\n")
for FILE in FILES:
  for iev,event in enumerate(FILE):
    NEvents.Fill(0)
    event.getByLabel(mu12Label, mu12)
    event.getByLabel(mu3Label, mu3)
    event.getByLabel(metsLabel, mets)
    event.getByLabel(tausLabel, taus)
    event.getByLabel(jetsLabel, jets)

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
     
    bestdR = 10000
    foundTau = False
    for i, tau in enumerate(taus.product() ):
      if tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") <= -0.5: continue
      for i, mu in enumerate(mu3.product() ):
        dr = DeltaR(tau.p4(), mu.p4() )
        #print "dr=", dr ,"\tbestdR=", bestdR
        if dr < bestdR:
          foundTau = True
          bestdR = dr
          tm = mu
          th = tau
          mva = tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT")
      bestTauJetdR = 100
      for i, jet in enumerate(jets.product() ):
        dr = DeltaR(tau.p4(), jet.p4() )
        if dr < bestTauJetdR:
          bestTauJetdR = dr
          tj = jet
    
    if foundTau == False: 
      continue
    NEvents.Fill(2)
    diTau = ROOT.TLorentzVector()
    diMu = ROOT.TLorentzVector()
    diMu = m1.p4() + m2.p4()
    diTau = tm.p4() + th.p4()
    met = mets.product().front()
    metp4 = KinematicFitter.MET('met', met.et(), met.phi(), met.getSignificanceMatrix()(0,0), met.getSignificanceMatrix()(0,1), met.getSignificanceMatrix()(1,0), met.getSignificanceMatrix()(1,1))

    ################
    # GET WEIGHTS  #
    ################
    # Pileup and Summed Weights
    nTrueVertices = -1
    if ISMC:
      event.getByLabel(pileupLabel, pileup)
      event.getByLabel(genEventLabel, genEvent)
      for cand in pileup.product():
        if cand.getBunchCrossing() == 0:
          nTrueVertices=cand.getTrueNumInteractions()
          break
      if nTrueVertices != -1:
        binx =  Pileup_.GetXaxis().FindBin(nTrueVertices)
        pileupWeight = Pileup_.GetBinContent(binx)
      if nTrueVertices > 35: nTrueVertices = 35 

      eventGenWeight = 1 #genEvent.weight()
      genWeight = eventGenWeight * LUMI * XSEC * .001 / SUMMEDWEIGHTS

      m1IDweight  = 1.0
      m1ISOweight = 1.0
      m2IDweight  = 1.0
      m2ISOweight = 1.0
      m3IDweight  = 1.0
      m1TrackWeight = 1.0
      m2TrackWeight = 1.0
      m3TrackWeight = 1.0
      LumiFraction_GH = 16.1 / 35.9
      LumiFraction_BF = 19.8 / 35.9 
    
      #mu1 SF
      binxIDs_BF = IDsWeight_BToF.GetXaxis().FindBin( m1.pt())
      binyIDs_BF = IDsWeight_BToF.GetYaxis().FindBin( math.fabs(m1.eta() ))
      binxIDs_GH = IDsWeight_GH.GetXaxis().FindBin(   m1.pt())
      binyIDs_GH = IDsWeight_GH.GetYaxis().FindBin(   math.fabs(m1.eta() ))
      binxISOs_BF = ISOsWeight_BToF.GetXaxis().FindBin( m1.pt())
      binyISOs_BF = ISOsWeight_BToF.GetYaxis().FindBin( math.fabs(m1.eta() ))
      binxISOs_GH = ISOsWeight_GH.GetXaxis().FindBin(   m1.pt())
      binyISOs_GH = ISOsWeight_GH.GetYaxis().FindBin(   math.fabs(m1.eta() ))
      if binxIDs_BF == 7:
        binxIDs_BF = 6
        binxIDs_GH = 6
        binxISOs_BF = 6
        binxISOs_GH = 6
      m1IDweight = IDsWeight_BToF.GetBinContent(  binxIDs_BF, binyIDs_BF) * LumiFraction_BF + LumiFraction_GH * IDsWeight_GH.GetBinContent(    binxIDs_GH, binyIDs_GH)    
      if DeltaR(m1.p4(), m2.p4()) > 0.4:
        m1ISOweight = ISOsWeight_BToF.GetBinContent(  binxISOs_BF, binyISOs_BF) * LumiFraction_BF + LumiFraction_GH * ISOsWeight_GH.GetBinContent(    binxISOs_GH, binyISOs_GH)

      # Track SF
      for i in TrackCorr_eta: 
        if math.fabs(m1.eta()) >= (i[0] - i[3]) and math.fabs(m1.eta()) <= (i[0] + i[2]):
          m1TrackWeight *= i[1]
          break
      for i in TrackCorr_vtx:
        if nTrueVertices >= (i[0] - i[3]) and nTrueVertices <= (i[0] + i[2]):
          m1TrackWeight *= i[1]
          break
  
      #mu2 SF
      if m2.pt() > 20.0:
        binxIDs_BF = IDsWeight_BToF.GetXaxis().FindBin( m2.pt() )
        binyIDs_BF = IDsWeight_BToF.GetYaxis().FindBin( math.fabs(m2.eta() ))
        binxIDs_GH = IDsWeight_GH.GetXaxis().FindBin(   m2.pt() )
        binyIDs_GH = IDsWeight_GH.GetYaxis().FindBin(   math.fabs(m2.eta() ))
        binxISOs_BF = ISOsWeight_BToF.GetXaxis().FindBin( m2.pt() )
        binyISOs_BF = ISOsWeight_BToF.GetYaxis().FindBin( math.fabs(m2.eta() ))
        binxISOs_GH = ISOsWeight_GH.GetXaxis().FindBin(   m2.pt() )
        binyISOs_GH = ISOsWeight_GH.GetYaxis().FindBin(   math.fabs(m2.eta() ))
        if binxIDs_BF == 7:
          binxIDs_BF = 6
          binxIDs_GH = 6
          binxISOs_BF = 6
          binxISOs_GH = 6
        m2IDweight = IDsWeight_BToF.GetBinContent(  binxIDs_BF, binyIDs_BF) * LumiFraction_BF + LumiFraction_GH * IDsWeight_GH.GetBinContent(    binxIDs_GH, binyIDs_GH)
        if DeltaR(m1.p4(), m2.p4()) > 0.4:
          m2ISOweight = ISOsWeight_BToF.GetBinContent(  binxISOs_BF, binyISOs_BF) * LumiFraction_BF + LumiFraction_GH * ISOsWeight_GH.GetBinContent(    binxISOs_GH, binyISOs_GH)
      else:
        binyIDs_BF  = IDsWeight_LowPt_ID.GetYaxis().FindBin( m2.pt())
        binxIDs_BF  = IDsWeight_LowPt_ID.GetXaxis().FindBin( math.fabs(m2.eta()))
        binyISOs_BF = ISOsWeight_LowPt_ISO.GetYaxis().FindBin( m2.pt())
        binxISOs_BF = ISOsWeight_LowPt_ISO.GetXaxis().FindBin( math.fabs(m2.eta()))
        m2IDweight  = IDsWeight_LowPt_ID.GetBinContent(  binxIDs_BF, binyIDs_BF)
        if DeltaR(m1.p4(), m2.p4()) > 0.4:
          m2ISOweight = ISOsWeight_LowPt_ISO.GetBinContent(  binxISOs_BF, binyISOs_BF)

      # Track SF
      for i in TrackCorr_eta:
        if math.fabs(m2.eta()) >= (i[0] - i[3]) and math.fabs(m2.eta()) <= (i[0] + i[2]):
          m2TrackWeight *= i[1]
          break
      for i in TrackCorr_vtx:
        if nTrueVertices >= (i[0] - i[3]) and nTrueVertices <= (i[0] + i[2]):
          m2TrackWeight *= i[1]
          break

      #mu3 SF
      if tm.pt() > 20.0:
        binxIDs_BF = IDsWeight_BToF.GetXaxis().FindBin( tm.pt() )
        binyIDs_BF = IDsWeight_BToF.GetYaxis().FindBin( math.fabs(tm.eta() ))
        binxIDs_GH = IDsWeight_GH.GetXaxis().FindBin(   tm.pt() )
        binyIDs_GH = IDsWeight_GH.GetYaxis().FindBin(   math.fabs(tm.eta() ))
        if binxIDs_BF == 7:
          binxIDs_BF = 6
          binxIDs_GH = 6
        m3IDweight = IDsWeight_BToF.GetBinContent(  binxIDs_BF, binyIDs_BF) * LumiFraction_BF + LumiFraction_GH * IDsWeight_GH.GetBinContent(    binxIDs_GH, binyIDs_GH)
      else:
        binyIDs_BF = IDsWeight_LowPt_ID.GetYaxis().FindBin( tm.pt())
        binxIDs_BF = IDsWeight_LowPt_ID.GetXaxis().FindBin( math.fabs(tm.eta()))
        m3IDweight = IDsWeight_LowPt_ID.GetBinContent(  binxIDs_BF, binyIDs_BF)

      # Track SF
      for i in TrackCorr_eta:
        if math.fabs(tm.eta()) >= (i[0] - i[3]) and math.fabs(tm.eta()) <= (i[0] + i[2]):
          m3TrackWeight *= i[1]
          break
      for i in TrackCorr_vtx:
        if nTrueVertices >= (i[0] - i[3]) and nTrueVertices <= (i[0] + i[2]):
          m3TrackWeight *= i[1]
          break
  
      # Trigger SF
      binxTrigger_BF = TriggerWeight_BToF.GetXaxis().FindBin( m1.pt())
      binyTrigger_BF = TriggerWeight_BToF.GetYaxis().FindBin( math.fabs(m1.eta()))
      binxTrigger_GH = TriggerWeight_GH.GetXaxis().FindBin(   m1.pt())
      binyTrigger_GH = TriggerWeight_GH.GetYaxis().FindBin(   math.fabs(m1.eta()))
      Trigger_weight = LumiFraction_BF *  TriggerWeight_BToF.GetBinContent( binxTrigger_BF, binyTrigger_BF) +  LumiFraction_GH *  TriggerWeight_GH.GetBinContent(   binxTrigger_GH, binyTrigger_GH)

      #Tau SF
      if th.pt() > 20.0: tauMedSF = 0.97
      else: tauMedSF = 0.70

      #BTag SF
      flavor = tj.hadronFlavour()
      if flavor == 5:   csvSF = reader.eval(0, tj.eta(), tj.pt() )
      elif flavor == 4: csvSF = reader.eval(1, tj.eta(), tj.pt() )
      else:             csvSF = reader.eval(2, tj.eta(), tj.pt() )
      if tj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > .8484: csvWeight = 1 - csvSF
      else: csvWeight = 1
      print 'tj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")= ', tj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), '\ncsvWeight=', csvWeight, '\tcsvSF=', csvSF

    elif APPLYFAKERATE:
      xaxis = histFR_.GetXaxis()
      yaxis = histFR_.GetYaxis()
      binx = xaxis.FindBin(th.pt()  )
      biny = yaxis.FindBin(math.fabs(th.eta()) )
      rate = histFR_.GetBinContent(binx, biny)
      print "rate= " , rate , "\tbinx=" ,  binx , "\tbiny= " , biny 
      while rate == 0 and binx > 0:
        print  "Bin(", binx, ",", biny, "_=0. New bin is (", binx-1, ",", biny , ")"
        binx -= 1
        rate = histFR_.GetBinContent(binx, biny)
        print " with new rate=" , rate 
      fakeRateWeight = rate / (1 - rate)
         
    #if ISMC: totalWeight = pileupWeight * genWeight * csvWeight * tauMedSF  * m1IDweight * m2IDweight * m3IDweight * m1ISOweight * m2ISOweight * Trigger_weight
    if ISMC: totalWeight = pileupWeight * genWeight * csvWeight * tauMedSF * m1IDweight * m2IDweight * m3IDweight * m1ISOweight * m2ISOweight * Trigger_weight
    elif APPLYFAKERATE: totalWeight = fakeRateWeight
    else: totalWeight = 1

    m1p4 = KinematicFitter.Muon(       'm1', m1.pt(), m1.eta(), m1.phi(), m1.energy())
    m2p4 = KinematicFitter.Muon(       'm2', m2.pt(), m2.eta(), m2.phi(), m2.energy())
    tmp4 = KinematicFitter.MuonTau(    'tm', tm.pt(), tm.eta(), tm.phi(), tm.energy())
    thp4 = KinematicFitter.HadronicTau('th', th.pt(), th.eta(), th.phi(), th.energy(), th.decayMode())
    m1p4.setErrors(0.01*m1p4.E(),0,0)
    m2p4.setErrors(0.01*m2p4.E(),0,0)
  
    kinfit = KinematicFitter.KinematicFitter()
    kinfit.addParticle('m1',m1p4)
    kinfit.addParticle('m2',m2p4)
    kinfit.addParticle('tm',tmp4)
    kinfit.addParticle('th',thp4)
    kinfit.addParticle('met',metp4)
    unc = 0.01
    if math.fabs(diMu.Eta() ) <= 1.2: unc = 0.018
    kinfit.addComposite('amm','m1','m2',uncertainty=unc*diMu.M())
    kinfit.addComposite('att','tm','th')
    kinfit.addComposite('h','m1','m2','tm','th')
    kinfit.addMassConstraint('amm','att',0)
    kinfit.setMinimizationParticles('th')
    kinfit.fit()     

    fourBody = m1.p4() + m2.p4() + tm.p4() + th.p4()
    hComp = kinfit.getComposite("h")
    #print type(hComp), hComp.M()
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
    TotalWeight.Fill( totalWeight) 
    hKF_ChiSq.Fill( kinfit.getChi2(), totalWeight )
    hKF_Mass.Fill( hComp.M(), totalWeight )

    x.setVal(  diMu.M() )
    y.setVal(  diTau.M() )
    y1.setVal( fourBody.M() )
    y2.setVal( hComp.M() )
    chi.setVal( kinfit.getChi2() )
    w.setVal( totalWeight )
    mumumass_dataset.add(         ROOT.RooArgSet(x,w) )
    mumutautaumass_dataset.add(   ROOT.RooArgSet(x,y,w) )
    mumufourBodymass_dataset.add( ROOT.RooArgSet(x,y1,w) )
    mumufourBodyKinFitmass_dataset.add(ROOT.RooArgSet(x,y2,chi,w) )

    diMudR = DeltaR(m1.p4(), m2.p4() )
    lumiWeight = LUMI * XSEC * .001 / SUMMEDWEIGHTS
    eventInfo = str(event.eventAuxiliary().run()) + "," + str(event.eventAuxiliary().luminosityBlock()) + "," + str(event.eventAuxiliary().event()) + ","
    m1Info  = str(m1.pt()) + "," + str(m1.eta()) + "," + str(m1.phi()) + "," + str(0) + "," + str(0) + ","
    m2Info  = str(m2.pt()) + "," + str(m2.eta()) + "," + str(m2.phi()) + "," + str(0) + "," + str(0) + ","
    tmInfo  = str(tm.pt()) + "," + str(tm.eta()) + "," + str(tm.phi()) + "," + str(0) + "," + str(0) + ","
    thInfo  = str(th.pt()) + "," + str(th.eta()) + "," + str(th.phi()) + "," + str(0) + "," + str(tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw")) + "," 
    thInfo2 = str(tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT")) + "," + str(tau.decayMode()) + ","
    aInfo = str(diMu.M()) + "," + str(diMudR) + "," + str(diTau.M()) + "," + str(bestdR) + ","
    hInfo = str(fourBody.M())  + "," + str(hComp.M()) + "," + str(kinfit.getChi2()) + ","     
    weightInfo  = str(eventGenWeight ) + "," + str(pileupWeight) + "," + str(Trigger_weight) + "," + str(m1IDweight) + "," + str(m1ISOweight) + "," + str(m2IDweight) + "," + str(m2ISOweight) + ","
    weightInfo2 = str(m3IDweight) + "," + str(m1TrackWeight) + "," + str(m2TrackWeight) + "," + str(m3TrackWeight) + "," + str(lumiWeight) + "," + str(csvWeight) + "," + str(totalWeight)
    csvFile.write(eventInfo + m1Info + m2Info + tmInfo + thInfo + thInfo2 + aInfo + hInfo + weightInfo + weightInfo2 + "\n")

outputFile.cd()
mumu_Mass.Write()
tautau_Mass.Write()
h_Mass.Write()
hKF_Mass.Write()
hKF_ChiSq.Write()
TotalWeight.Write()
DiTauMassVSMVA.Write()
DiTauMassVSTauMu3dR.Write()
DiTauMassVSDiMudR.Write()
DiTauMassVSTauHMu1dR.Write()
DiTauMassVSmu3Pt.Write()
DiTauMassVStauPt.Write()
DiTauMassVSDiTauPt.Write()
DiTauMassVSDiMuMass.Write()
mumumass_dataset.Write()
mumutautaumass_dataset.Write()
mumufourBodymass_dataset.Write()
mumufourBodyKinFitmass_dataset.Write()
outputFile.Write()
outputFile.Close()
csvFile.close()
