import sys
import os
import json

from math import floor

import ROOT


class PileupWeights(object):

    def __init__(self,version):
        self.version = version
        path = '{0}/src/DevTools/Analyzer/data/pileup_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12.root'.format(os.environ['CMSSW_BASE'])
        self.scale = {}
        self.scale_up = {}
        self.scale_down = {}
        self.alt_scales = {}
        rootfile = ROOT.TFile(path)
        hist_scale = rootfile.Get('pileup_scale')
        for b in range(hist_scale.GetNbinsX()):
            self.scale[b] = hist_scale.GetBinContent(b+1)
        hist_scale = rootfile.Get('pileup_scale_up')
        for b in range(hist_scale.GetNbinsX()):
            self.scale_up[b] = hist_scale.GetBinContent(b+1)
        hist_scale = rootfile.Get('pileup_scale_down')
        for b in range(hist_scale.GetNbinsX()):
            self.scale_down[b] = hist_scale.GetBinContent(b+1)
        for xsec in [65000,66000,67000,68000,69000,70000,71000,72000,73000,74000,75000,76000,77000,78000,79000,80000]:
            self.alt_scales[xsec] = {}
            hist_scale = rootfile.Get('pileup_scale_{0}'.format(xsec))
            for b in range(hist_scale.GetNbinsX()):
                self.alt_scales[xsec][b] = hist_scale.GetBinContent(b+1)
        rootfile.Close()

    def alt_weight(self, rtrow,xsec):
        if rtrow.nTrueVertices < 0:
            return 1
        else:
            if xsec in self.alt_scales:
                val = self.alt_scales[xsec][int(floor(rtrow.nTrueVertices))]
            else:
                val = 1
            return val

    def weight(self, rtrow):
        if rtrow.nTrueVertices < 0:
            return [1,1,1]
        else:
            val = self.scale[int(floor(rtrow.nTrueVertices))]
            up = self.scale_up[int(floor(rtrow.nTrueVertices))]
            down = self.scale_down[int(floor(rtrow.nTrueVertices))]
            return [val,up,down]
