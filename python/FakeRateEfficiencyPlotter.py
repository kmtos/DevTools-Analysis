import string
import numpy as np
import csv
import pandas as pd
import sys
from collections import defaultdict
from ROOT import TCanvas, TGraphErrors, TLegend, TFile
from ROOT import gROOT
import ROOT
from math import sin
from array import array
import CMS_lumi as CMS_lumi

c1 = TCanvas( 'c1', 'A Simple Graph Example', 800,600)
c1.SetTopMargin(0.06);
c1.SetRightMargin(0.03);
c1.SetBottomMargin(0.12);
c1.SetLeftMargin(0.11);

efficiencies = { 'Cleaned':     {'125': {}, '300': {}, '750': {} },
                 'NotCleaned': {'125': {}, '300': {}, '750': {} }
               }

with open("/home/kyletos/Downloads/RootFiles/ABCD_Tests/FakeRateEfficiency/RatioValues.out", "r") as f:  
  for line in f:
    hmass, amass, region, eff, count = line.split()
    efficiencies[region][str(hmass)][str(amass)] = (eff, count)

#for k,v in efficiencies.items():
#  for k1,v1, in efficiencies[k].items():
#    for k2,v2, in efficiencies[k][k1].items():
#      print k, k1, k2, v2

temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["Cleaned"]["125"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1]) ) )
temp.sort(key=lambda tup: tup[0])
for i in temp: 
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): 
  print erry[i] , yVals[i] * (1 - yVals[i]), (yVals[i] * (1 - yVals[i]) / erry[i]), (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
  erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g125C = TGraphErrors(9,xVals,yVals,errx,erry)
g125C.SetLineColor(ROOT.kGreen+2)
g125C.SetLineWidth(3)
g125C.SetMarkerColor(ROOT.kGreen+2)
g125C.SetLineStyle(1)
g125C.SetMarkerSize(1)
g125C.SetMarkerStyle(8)
g125C.SetTitle("")
g125C.GetXaxis().SetTitle("m_{a} GeV")
#g125C.GetXaxis().SetTitleSize(1)
g125C.GetYaxis().SetTitle("Reconstruction Efficiency")
#g125C.GetYaxis().SetTitleSize(1)

temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["Cleaned"]["300"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1])) )
temp.sort(key=lambda tup: tup[0])
for i in temp:
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g300C = TGraphErrors(9,xVals,yVals,errx,erry)
g300C.SetLineColor(ROOT.kBlue)
g300C.SetLineWidth(3)
g300C.SetMarkerColor(ROOT.kBlue)
g300C.SetLineStyle(1)
g300C.SetMarkerSize(1)
g300C.SetMarkerStyle(8)
g300C.SetTitle("")
g300C.GetXaxis().SetTitle("m_{a} GeV")
#g300C.GetXaxis().SetTitleSize(1)
g300C.GetYaxis().SetTitle("Reconstruction Efficiency")
#g300C.GetYaxis().SetTitleSize(1)

temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["Cleaned"]["750"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1])) )
temp.sort(key=lambda tup: tup[0])
for i in temp:
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g750C = TGraphErrors(9,xVals,yVals,errx,erry)
g750C.SetLineColor(ROOT.kRed+1)
g750C.SetLineWidth(3)
g750C.SetMarkerColor(ROOT.kRed+1)
g750C.SetLineStyle(1)
g750C.SetMarkerSize(1)
g750C.SetMarkerStyle(8)
g750C.SetTitle("")
g750C.GetXaxis().SetTitle("m_{a} GeV")
#g750C.GetXaxis().SetTitleSize(1)
g750C.GetYaxis().SetTitle("Reconstruction Efficiency")
#g750C.GetYaxis().SetTitleSize(1)


temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["NotCleaned"]["125"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1])) )
temp.sort(key=lambda tup: tup[0])
for i in temp:
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g125N = TGraphErrors(9,xVals,yVals,errx,erry)
g125N.SetLineColor(ROOT.kGreen+2)
g125N.SetLineWidth(3)
g125N.SetMarkerColor(ROOT.kGreen+2)
g125N.SetLineStyle(2)
g125N.SetMarkerSize(1)
g125N.SetMarkerStyle(8)
g125N.SetTitle("")
g125N.GetXaxis().SetTitle("m_{a} GeV")
#g125N.GetXaxis().SetTitleSize(1)
g125N.GetYaxis().SetTitle("Reconstruction Efficiency")
#g125N.GetYaxis().SetTitleSize(1)

temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["NotCleaned"]["300"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1])) )
temp.sort(key=lambda tup: tup[0])
for i in temp:
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g300N = TGraphErrors(9,xVals,yVals,errx,erry)
g300N.SetLineColor(ROOT.kBlue)
g300N.SetLineWidth(3)
g300N.SetMarkerColor(ROOT.kBlue)
g300N.SetLineStyle(2)
g300N.SetMarkerSize(1)
g300N.SetMarkerStyle(8)
g300N.SetTitle("")
g300N.GetXaxis().SetTitle("m_{a} GeV")
#g300N.GetXaxis().SetTitleSize(1)
g300N.GetYaxis().SetTitle("Reconstruction Efficiency")
#g300N.GetYaxis().SetTitleSize(1)

temp = []
xVals = array( 'd' )
yVals = array( 'd' )
errx = array( 'd' )
erry = array( 'd' )
for k2,v2, in efficiencies["NotCleaned"]["750"].items():
  temp.append( (int(k2),float(v2[0]), float(v2[1])) )
temp.sort(key=lambda tup: tup[0])
for i in temp:
  xVals.append(i[0])
  yVals.append(i[1])
  errx.append(0)
  erry.append(i[2])
for i in range(len(erry) ): erry[i] = (yVals[i] * (1 - yVals[i]) / erry[i]) **(1/2.0)
g750N = TGraphErrors(9,xVals,yVals,errx,erry)
g750N.SetLineColor(ROOT.kRed+1)
g750N.SetLineWidth(3)
g750N.SetMarkerColor(ROOT.kRed+1)
g750N.SetLineStyle(2)
g750N.SetMarkerSize(1)
g750N.SetMarkerStyle(8)
g750N.SetTitle("")
g750N.GetXaxis().SetTitle("m_{a} GeV")
#g750N.GetXaxis().SetTitleSize(1.1)
g750N.GetYaxis().SetTitle("Reconstruction Efficiency")
#g750N.GetYaxis().SetTitleSize(1.1)

leg = TLegend(.3,.18,.55,.4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.SetHeader("Muon Cleaned HPS")
leg.AddEntry(g125C,"m_{h} = 125 GeV","L")
leg.AddEntry(g300C,"m_{h} = 300 GeV","L")
leg.AddEntry(g750C,"m_{h} = 750 GeV","L")

leg2 = TLegend(.6,.18,.85,.4)
leg2.SetBorderSize(0)
leg2.SetFillColor(0)
leg2.SetFillStyle(0)
leg2.SetTextFont(42)
leg2.SetTextSize(0.035)
leg2.SetHeader("HPS")
leg2.AddEntry(g125N,"m_{h} = 125 GeV","L")
leg2.AddEntry(g300N,"m_{h} = 300 GeV","L")
leg2.AddEntry(g750N,"m_{h} = 750 GeV","L")

mg = ROOT.TMultiGraph()
mg.Add(g125C)
mg.Add(g300C)
mg.Add(g750C)
mg.Add(g125N)
mg.Add(g300N)
mg.Add(g750N)

c1.cd()
#g125C.Draw("ALP")
#g300C.Draw("LP")
#g750C.Draw("LP")
#g125N.Draw("LP")
#g300N.Draw("LP")
#g750N.Draw("LP")
mg.Draw("ALP")
mg.GetXaxis().SetTitle("m_{a} GeV");
mg.GetXaxis().SetTitleSize(0.05);
mg.GetYaxis().SetTitle("Reconstruction Efficiency");
mg.GetYaxis().SetTitleSize(0.05);
mg.GetYaxis().SetRangeUser(0.,1.0);
c1.SetTickx()
c1.SetTicky()
leg.Draw("SAME")
leg2.Draw("SAME")

CMS_lumi.cmsText = 'CMS' 
CMS_lumi.writeExtraText = True 
CMS_lumi.extraText = "Preliminary Simulation" 
CMS_lumi.lumi_13TeV = ''
CMS_lumi.cmsTextSize =0.4 
CMS_lumi.CMS_lumi(c1,123,11)
f = ROOT.TFile("/home/kyletos/Downloads/RootFiles/ABCD_Tests/FakeRateEfficiency/CleanedEfficiencies.root", "recreate")
f.cd()
c1.Write()
f.Write()
f.Close()
c1.SaveAs("/home/kyletos/Downloads/RootFiles/ABCD_Tests/FakeRateEfficiency/CleanedEfficiencies.pdf", "pdf")
c1.Print("/home/kyletos/Downloads/RootFiles/ABCD_Tests/FakeRateEfficiency/CleanedEfficiencies.png")
