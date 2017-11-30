import sys
import os
import json
import logging
from math import floor

import ROOT


class PileupWeights(object):

    def __init__(self,version,profile=None):
        # TODO: support mixed data/mc campaigns
        self.version = version
        if version == '76X':
            if not profile: profile = 'PU25nsData2015v1'
        else:
            if not profile: profile = 'PUMoriond17'

        paths = {
            ('76X','PU25nsData2015v1') : 'pileup_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12.root',
            ('80X','PUMoriond17')      : 'pileup_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6.root',
        }

        self.skip = (version,profile) not in paths
        if self.skip:
            logging.warning('No pileup found for {0} {1}. No weights will be applied.'.format(version,profile))
            return

        path = '{0}/src/DevTools/Analyzer/data/{1}'.format(os.environ['CMSSW_BASE'],paths[(version,profile)])

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
        for xsec in [60000,61000,62000,63000,64000,65000,66000,67000,68000,69000,70000,71000,72000,73000,74000,75000,76000,77000,78000,79000,80000]:
            self.alt_scales[xsec] = {}
            hist_scale = rootfile.Get('pileup_scale_{0}'.format(xsec))
            for b in range(hist_scale.GetNbinsX()):
                self.alt_scales[xsec][b] = hist_scale.GetBinContent(b+1)
        rootfile.Close()

    def alt_weight(self,event,xsec):
        if self.skip: return 1
        vert = event.nTrueVertices()
        isData = event.isData()>0.5
        if vert < 0 or isData:
            return 1
        else:
            if xsec in self.alt_scales:
                val = self.alt_scales[xsec][int(floor(vert))]
            else:
                val = 1
            return val

    def weight(self, event):
        if self.skip: return [1,1,1]
        vert = event.nTrueVertices()
        isData = event.isData()>0.5
        if vert < 0 or isData:
            return [1,1,1]
        else:
            val = self.scale[int(floor(vert))]
            up = self.scale_up[int(floor(vert))]
            down = self.scale_down[int(floor(vert))]
            return [val,up,down]
