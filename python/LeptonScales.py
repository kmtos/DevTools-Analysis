import os
import sys
import logging

import ROOT

from DevTools.Utilities.utilities import *

class LeptonScales(object):
    '''Class to access the lepton scales for a given ID.'''

    def __init__(self,version):
        self.version = version

        # EGamma POG
        self.egamma_pog_scales = {}
        self.egamma_pog_rootfiles = {}
        for cbid in ['Veto','Loose','Medium','Tight']:
            name = 'Cutbased{0}'.format(cbid)
            path = '{0}/src/DevTools/Analyzer/data/CutBasedID_{1}WP_76X_18Feb.txt_SF2D.root'.format(os.environ['CMSSW_BASE'],cbid)
            if self.version=='80X':
                path = '{0}/src/DevTools/Analyzer/data/egammaEffi.txt_SF2D_CBID{1}_2016.root'.format(os.environ['CMSSW_BASE'],cbid)
            self.egamma_pog_rootfiles[name] = ROOT.TFile(path)
            self.egamma_pog_scales[name] = self.egamma_pog_rootfiles[name].Get('EGamma_SF2D')

        if self.version=='80X':
            name = 'GSFTracking'
            path = '{0}/src/DevTools/Analyzer/data/gsf_tracking_scalefactor_2016.root'.format(os.environ['CMSSW_BASE'])
            self.egamma_pog_rootfiles[name] = ROOT.TFile(path)
            self.egamma_pog_scales[name] = self.egamma_pog_rootfiles[name].Get('EGamma_SF2D')

        # Muon POG
        self.muon_pog_scales = {}
        idpath = '{0}/src/DevTools/Analyzer/data/MuonID_Z_RunCD_Reco76X_Feb15.root'.format(os.environ['CMSSW_BASE'])
        isopath = '{0}/src/DevTools/Analyzer/data/MuonIso_Z_RunCD_Reco76X_Feb15.root'.format(os.environ['CMSSW_BASE'])
        if self.version=='80X':
            idpath = '{0}/src/DevTools/Analyzer/data/MuonID_Z_RunBCD_prompt80X_7p65.root'.format(os.environ['CMSSW_BASE'])
            isopath = '{0}/src/DevTools/Analyzer/data/MuonIso_Z_RunBCD_prompt80X_7p65.root'.format(os.environ['CMSSW_BASE'])
        self.muon_pog_id_rootfile = ROOT.TFile(idpath)
        self.muon_pog_iso_rootfile = ROOT.TFile(isopath)
        for mid in ['LooseID','MediumID','SoftID','TightID']:
            idname = mid + 'andIPCut' if mid=='TightID' else mid
            self.muon_pog_scales[mid] = self.muon_pog_id_rootfile.Get('MC_NUM_{0}_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio'.format(idname))
            for iso in ['LooseRelIso','TightRelIso']:
                if mid=='SoftID': continue
                if mid=='LooseID' and iso=='LooseRelIso': continue
                name = '{0}{1}'.format(iso,mid)
                if self.version=='80X': mid = 'TightID' # others not included yet
                self.muon_pog_scales[name] = self.muon_pog_iso_rootfile.Get('MC_NUM_{0}_DEN_{1}_PAR_pt_spliteta_bin1/abseta_pt_ratio'.format(iso,mid))


        # 80X HIP efficiency
        # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Tracking_efficiency_provided_by
        self.tracking_pog_scales = {}
        if self.version=='80X':
            scalepathless10 = '{0}/src/DevTools/Analyzer/data/muon_tracking_eff_2016_lesspt10.root'.format(os.environ['CMSSW_BASE'])
            scalepathgreater10 = '{0}/src/DevTools/Analyzer/data/muon_tracking_eff_2016_greaterpt10.root'.format(os.environ['CMSSW_BASE'])
            rootless = ROOT.TFile(scalepathless10)
            rootgreater = ROOT.TFile(scalepathgreater10)
            self.tracking_pog_scales['TrackingLessPt10NVtx']    = self.__parseAsymmErrors(rootless.Get('ratio_vtx'))
            self.tracking_pog_scales['TrackingLessPt10Eta']     = self.__parseAsymmErrors(rootless.Get('ratio_eta'))
            self.tracking_pog_scales['TrackingGreaterPt10NVtx'] = self.__parseAsymmErrors(rootgreater.Get('ratio_vtx'))
            self.tracking_pog_scales['TrackingGreaterPt10Eta']  = self.__parseAsymmErrors(rootgreater.Get('ratio_eta'))
            rootless.Close()
            rootgreater.Close()
        

        # HZZ
        self.muon_hzz_scales = {}
        idpath = '{0}/src/DevTools/Analyzer/data/HZZ_muon_76X.root'.format(os.environ['CMSSW_BASE'])
        self.muon_hzz_rootfile = ROOT.TFile(idpath)
        self.muon_hzz_scales['HZZTight'] = self.muon_hzz_rootfile.Get('FINAL')

        self.electron_hzz_scales = {}
        idpath = '{0}/src/DevTools/Analyzer/data/HZZ_electron_76X.root'.format(os.environ['CMSSW_BASE'])
        self.electron_hzz_rootfile = ROOT.TFile(idpath)
        self.electron_hzz_scales['HZZTight'] = self.electron_hzz_rootfile.Get('hScaleFactors_IdIsoSip')

        # private 80X
        
        # electron
        self.private_electron_80X = {}
        path = '{0}/src/DevTools/Analyzer/data/scalefactors_electron_2016.root'.format(os.environ['CMSSW_BASE'])
        self.private_electron_80X_rootfile = ROOT.TFile(path)
        for idName in ['CutBasedIDVeto','CutBasedIDLoose','CutBasedIDMedium','CutBasedIDTight']:
            self.private_electron_80X[idName] = self.private_electron_80X_rootfile.Get(idName)

        # muon
        self.private_muon_80X = {}
        path = '{0}/src/DevTools/Analyzer/data/scalefactors_muon_2016.root'.format(os.environ['CMSSW_BASE'])
        self.private_muon_80X_rootfile = ROOT.TFile(path)
        for idName in ['LooseID','LooseIsoFromLooseID','MediumID','LooseIsoFromMediumID','TightIsoFromMediumID','MediumIDICHEP','LooseIsoFromMediumIDICHEP','TightIsoFromMediumIDICHEP','TightID','TightIsoFromTightID']:
            self.private_muon_80X[idName] = self.private_muon_80X_rootfile.Get(idName)

    def __exit__(self, type, value, traceback):
        self.finish()

    def __del__(self):
        self.finish()

    def finish(self):
        for name,rootfile in self.egamma_pog_rootfiles.iteritems():
            rootfile.Close()
        self.muon_pog_id_rootfile.Close()
        self.muon_pog_iso_rootfile.Close()
        self.muon_hzz_rootfile.Close()
        self.electron_hzz_rootfile.Close()
        self.private_electron_80X_rootfile.Close()
        self.private_muon_80X_rootfile.Close()


    def __parseAsymmErrors(self,graph):
        vals = []
        x,y = ROOT.Double(0), ROOT.Double(0)
        for i in range(graph.GetN()):
            graph.GetPoint(i+1,x,y)
            val = {
                'x'        : float(x),
                'y'        : float(y),
                'errx_up'  : float(graph.GetErrorXhigh(i+1)),
                'errx_down': float(graph.GetErrorXlow(i+1)),
                'erry_up'  : float(graph.GetErrorYhigh(i+1)),
                'erry_down': float(graph.GetErrorYlow(i+1)),
            }
            vals += [val]
        return vals
            

    def __getElectronScale(self,leptonId,cand):
        pt  = cand.pt()
        eta = cand.eta()
        sceta = cand.superClusterEta()
        if leptonId in self.egamma_pog_scales: # default to POG
            if pt>200: pt = 199.
            if pt<10: pt = 11.
            if 'Trig' in leptonId and pt<15: pt = 16.
            if 'GSFTracking' in leptonId and pt<20: pt = 21
            hist = self.egamma_pog_scales[leptonId]
            b = hist.FindBin(abs(sceta),pt) if self.version=='76X' else hist.FindBin(sceta,pt)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        elif self.version=='76X' and leptonId in self.electron_hzz_scales:
            if pt>200: pt = 199
            if pt<7: pt = 8
            hist = self.electron_hzz_scales[leptonId]
            b = hist.FindBin(pt,eta)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        elif self.version=='80X' and leptonId in self.private_electron_80X: # fall back to private
            if pt>200: pt = 199.
            if pt<10: pt = 11.
            hist = self.private_electron_80X[leptonId]
            b = hist.FindBin(pt,eta)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        else:
            val = 1.
            err = 0.
        return val, err

    def __getMuonScale(self,leptonId,cand):
        pt  = cand.pt()
        eta = cand.eta()
        if leptonId in self.muon_pog_scales: # default to pog
            if pt>120: pt = 119.
            if pt<20: pt = 21.
            hist = self.muon_pog_scales[leptonId]
            b = hist.FindBin(abs(eta),pt)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        elif self.version=='76X' and leptonId in self.muon_hzz_scales:
            if pt>80: pt = 79
            if pt<5: pt = 6
            hist = self.muon_hzz_scales[leptonId]
            b = hist.FindBin(pt,eta)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        elif self.version=='80X' and leptonId in self.private_muon_80X: # fallback to private
            if pt>200: pt = 199.
            if pt<10: pt = 11.
            hist = self.private_muon_80X[leptonId]
            b = hist.FindBin(pt,eta)
            val = hist.GetBinContent(b)
            err = hist.GetBinError(b)
        else:
            val = 1.
            err = 0.
        return val, err

    def __getTauScale(self,leptonId,cand):
        #pt  = cand.pt()
        #eta = cand.eta()
        return 1., 0. # simple recommendation, 6% error

    def __getMuonTrackingScale(self,cand):
        pt  = cand.pt()
        eta = cand.eta()
        nvtx = cand.tree.nTrueVertices
        if nvtx>35: nvtx = 35
        etaName = 'Tracking{0}Pt10Eta'.format('Less' if pt<10 else 'Greater')
        nvtxName = 'Tracking{0}Pt10NVtx'.format('Less' if pt<10 else 'Greater')
        valEta, errEta = (1., 0.)
        valNVtx, errNVtx = (1., 0.)
        if etaName in self.tracking_pog_scales:
            for valDict in self.tracking_pog_scales[etaName]:
                etaLow = valDict['x'] - valDict['errx_down']
                etaHigh = valDict['x'] + valDict['errx_up']
                if eta>=etaLow and eta<=etaHigh:
                    valEta = valDict['y']
                    errEta = (valDict['erry_up']+valDict['erry_down'])/2.
        if nvtxName in self.tracking_pog_scales:
            for valDict in self.tracking_pog_scales[nvtxName]:
                nvtxLow = valDict['x'] - valDict['errx_down']
                nvtxHigh = valDict['x'] + valDict['errx_up']
                if nvtx>=nvtxLow and nvtx<=nvtxHigh:
                    valNVtx = valDict['y']
                    errNVtx = (valDict['erry_up']+valDict['erry_down'])/2.
        val, err = prodWithError((valEta,errEta),(valNVtx,errNVtx))
        return val,err


    def getScale(self,leptonId,cand,doError=False):
        '''Get the scale to apply to MC (eff_data/eff_mc)'''
        if cand.collName=='electrons':
            idval = self.__getElectronScale(leptonId,cand)
            trackval = self.__getElectronScale('GSFTracking',cand)
            if cand.pt()< 20 and self.version=='80X':
                # additional 3% uncertainty
                # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
                val, err = prodWithError(idval,trackval,(1.,0.03))
            else:
                val, err = prodWithError(idval,trackval)
        elif cand.collName=='muons':
            if leptonId == 'TightIDTightIso':
                idname, isoname = ('TightID', 'TightRelIsoTightID')# if self.version=='76X' else ('TightID', 'TightIsoFromTightID')
            elif leptonId == 'MediumIDTightIso':
                idname, isoname = ('MediumID', 'TightRelIsoMediumID')# if self.version=='76X' else ('MediumIDICHEP', 'TightIsoFromMediumIDICHEP')
            elif leptonId == 'MediumIDLooseIso':
                idname, isoname = ('MediumID', 'LooseRelIsoMediumID')# if self.version=='76X' else ('MediumIDICHEP', 'LooseIsoFromMediumIDICHEP')
            else:
                idname, isoname = '', ''
            if idname and isoname:
                idval = self.__getMuonScale(idname,cand)
                isoval = self.__getMuonScale(isoname,cand)
                val, err = prodWithError(idval,isoval)
            else:
                val, err = self.__getMuonScale(leptonId,cand)
            valTrack, errTrack = self.__getMuonTrackingScale(cand)
            #val, err = prodWithError((val,err),(valTrack,errTrack))
            val, err = prodWithError((val,err),(valTrack,0. if errTrack!=errTrack else errTrack)) # bug with error from tgraphasymm, check for NaN
        elif cand.collName=='taus':
            val, err = self.__getTauScale(leptonId,cand)
        else:
            val, err = 1., 0.
        return (val, val+err, max([0.,val-err])) if doError else val
