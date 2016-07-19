import os
import sys
import logging

import ROOT

class LeptonScales(object):
    '''Class to access the lepton scales for a given ID.'''

    def __init__(self,version):
        self.version = version
        # EGamma POG
        self.egamma_pog_scales = {}
        self.egamma_pog_rootfiles = {}
        for cbid in ['Veto','Loose','Medium','Tight']:
            path = '{0}/src/DevTools/Analyzer/data/CutBasedID_{1}WP_76X_18Feb.txt_SF2D.root'.format(os.environ['CMSSW_BASE'],cbid)
            name = 'Cutbased{0}'.format(cbid)
            self.egamma_pog_rootfiles[name] = ROOT.TFile(path)
            self.egamma_pog_scales[name] = self.egamma_pog_rootfiles[name].Get('EGamma_SF2D')
        for wp in ['TrigWP80','TrigWP90']:
            path = '{0}/src/DevTools/Analyzer/data/ScaleFactor_GsfElectronToRECO_passing{1}.txt.egamma_SF2D.root'.format(os.environ['CMSSW_BASE'],wp)
            name = wp
            self.egamma_pog_rootfiles[name] = ROOT.TFile(path)
            self.egamma_pog_scales[name] = self.egamma_pog_rootfiles[name].Get('EGamma_SF2D')
        # Muon POG
        self.muon_pog_scales = {}
        idpath = '{0}/src/DevTools/Analyzer/data/MuonID_Z_RunCD_Reco76X_Feb15.root'.format(os.environ['CMSSW_BASE'])
        isopath = '{0}/src/DevTools/Analyzer/data/MuonIso_Z_RunCD_Reco76X_Feb15.root'.format(os.environ['CMSSW_BASE'])
        self.muon_pog_id_rootfile = ROOT.TFile(idpath)
        self.muon_pog_iso_rootfile = ROOT.TFile(isopath)
        for mid in ['LooseID','MediumID','SoftID','TightID']:
            idname = mid + 'andIPCut' if mid=='TightID' else mid
            self.muon_pog_scales[mid] = self.muon_pog_id_rootfile.Get('MC_NUM_{0}_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio'.format(idname))
            for iso in ['LooseRelIso','TightRelIso']:
                if mid=='SoftID': continue
                if mid=='LooseID' and iso=='LooseRelIso': continue
                name = '{0}{1}'.format(iso,mid)
                self.muon_pog_scales[name] = self.muon_pog_iso_rootfile.Get('MC_NUM_{0}_DEN_{1}_PAR_pt_spliteta_bin1/abseta_pt_ratio'.format(iso,mid))
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

    def __getElectronScale(self,leptonId,cand):
        pt  = cand.pt()
        eta = cand.eta()
        if self.version=='76X':
            if leptonId in self.egamma_pog_scales:
                if pt>200: pt = 199.
                if pt<10: pt = 11.
                if 'Trig' in leptonId and pt<15: pt = 16.
                hist = self.egamma_pog_scales[leptonId]
                b = hist.FindBin(abs(eta),pt)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            elif leptonId in self.electron_hzz_scales:
                if pt>200: pt = 199
                if pt<7: pt = 8
                hist = self.electron_hzz_scales[leptonId]
                b = hist.FindBin(pt,eta)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            else:
                val = 1.
        elif self.version=='80X':
            if leptonId in self.private_electron_80X:
                if pt>200: pt = 199.
                if pt<10: pt = 11.
                hist = self.private_electron_80X[leptonId]
                b = hist.FindBin(pt,eta)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            else:
                val = 1.
                err = 0.
        else:
            val = 1.
            err = 0.
        return val, err

    def __getMuonScale(self,leptonId,cand):
        pt  = cand.pt()
        eta = cand.eta()
        if self.version=='76X':
            if leptonId in self.muon_pog_scales:
                if pt>120: pt = 119.
                if pt<20: pt = 21.
                hist = self.muon_pog_scales[leptonId]
                b = hist.FindBin(abs(eta),pt)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            elif leptonId in self.muon_hzz_scales:
                if pt>80: pt = 79
                if pt<5: pt = 6
                hist = self.muon_hzz_scales[leptonId]
                b = hist.FindBin(pt,eta)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            else:
                val = 1.
                err = 0.
        elif self.version=='80X':
            if leptonId in self.private_muon_80X:
                if pt>200: pt = 199.
                if pt<10: pt = 11.
                hist = self.private_muon_80X[leptonId]
                b = hist.FindBin(pt,eta)
                val = hist.GetBinContent(b)
                err = hist.GetBinError(b)
            else:
                val = 1.
                err = 0.
        else:
            val = 1.
            err = 0.
        return val, err

    def __getTauScale(self,leptonId,cand):
        #pt  = cand.pt()
        #eta = cand.eta()
        return 1., 0. # simple recommendation, 6% error

    def getScale(self,leptonId,cand,doError=False):
        '''Get the scale to apply to MC (eff_data/eff_mc)'''
        if cand.collName=='electrons':
            val, err = self.__getElectronScale(leptonId,cand)
        elif cand.collName=='muons':
            if leptonId == 'TightIDTightIso':
                idname, isoname = ('TightID', 'TightRelIsoTightID') if self.version=='76X' else ('TightID', 'TightIsoFromTightID')
            elif leptonId == 'MediumIDTightIso':
                idname, isoname = ('MediumID', 'TightRelIsoMediumID') if self.version=='76X' else ('MediumIDICHEP', 'TightIsoFromMediumIDICHEP')
            elif leptonId == 'MediumIDLooseIso':
                idname, isoname = ('MediumID', 'LooseRelIsoMediumID') if self.version=='76X' else ('MediumIDICHEP', 'LooseIsoFromMediumIDICHEP')
            else:
                idname, isoname = '', ''
            if idname and isoname:
                idval = self.__getMuonScale(idname,cand)
                isoval = self.__getMuonScale(isoname,cand)
                val, err = idval*isoval
            else:
                val, err = self.__getMuonScale(leptonId,cand)
        elif cand.collName=='taus':
            val, err = self.__getTauScale(leptonId,cand)
        else:
            val, err = 1., 0.
        return (val, val+err, max([0.,val-err])) if doError else val
