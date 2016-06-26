import os
import sys
import logging
import math

import ROOT

import operator

def product(iterable):
    if not iterable: return 0. # case of empty list
    return reduce(operator.mul, iterable, 1)

class TriggerScales(object):
    '''Class to access the trigger scales for a given trigger configuration.'''

    def __init__(self,version):
        self.version = version
        ####################
        ### POG APPROVED ###
        ####################
        # single muon
        # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2
        singleMu_path = '{0}/src/DevTools/Analyzer/data/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root'.format(os.environ['CMSSW_BASE'])
        self.singleMu_rootfile = ROOT.TFile(singleMu_path)
        self.singleMu_efficiencies = {}
        for name in ['runC_IsoMu20_OR_IsoTkMu20', 'runC_Mu45_eta2p1', 'runC_Mu50',
                     'runD_IsoMu20_OR_IsoTkMu20_HLTv4p2', 'runD_IsoMu20_OR_IsoTkMu20_HLTv4p3',
                     'runD_Mu45_eta2p1', 'runD_Mu50']:
                directory = '{0}_PtEtaBins'.format(name)
                self.singleMu_efficiencies[name] = {
                    'MC'   : self.singleMu_rootfile.Get('{0}/efficienciesMC/pt_abseta_MC'.format(directory)),
                    'DATA' : self.singleMu_rootfile.Get('{0}/efficienciesDATA/pt_abseta_DATA'.format(directory)),
                }

        # double tau fits
        self.doubleTau_efficiencies = {
            'MC': {
                'val'    : {'m0' : 3.60274e01, 'sigma' : 5.89434e00, 'alpha' : 5.82870e00, 'n' : 1.83737e00, 'norm' : 9.58000e-01,},
                'errUp'  : {'m0' : 3.56012e01, 'sigma' : 5.97209e00, 'alpha' : 6.09604e00, 'n' : 1.68740e00, 'norm' : 9.87653e-01,},
                'errDown': {'m0' : 3.62436e01, 'sigma' : 5.58461e00, 'alpha' : 5.12924e00, 'n' : 2.05921e00, 'norm' : 9.32305e-01,},
            },
            'DATA': {
                'val'    : {'m0' : 3.45412e01, 'sigma' : 5.63353e00, 'alpha' : 2.49242e00, 'n' : 3.35896e00, 'norm' : 1.00000e00,},
                'errUp'  : {'m0' : 3.31713e01, 'sigma' : 5.66551e00, 'alpha' : 1.87175e00, 'n' : 8.07790e00, 'norm' : 1.00000e00,},
                'errDown': {'m0' : 3.56264e01, 'sigma' : 5.30711e00, 'alpha' : 2.81591e00, 'n' : 2.40649e00, 'norm' : 9.99958e-01,},
            },
        }

        ##############
        ### OTHERS ###
        ##############
        # HWW measurements: single mu, single e, double mu per leg, double e per leg
        # https://twiki.cern.ch/twiki/bin/view/CMS/HWW2015TriggerResults
        singleMu_path     = '{0}/src/DevTools/Analyzer/data/HWW_HLT_IsoMu20orIsoTkMu20_76X.txt'.format(os.environ['CMSSW_BASE'])
        doubleMuLeg1_path = '{0}/src/DevTools/Analyzer/data/HWW_HLT_Mu17_Mu8Leg1_76X.txt'.format(os.environ['CMSSW_BASE'])
        doubleMuLeg2_path = '{0}/src/DevTools/Analyzer/data/HWW_HLT_Mu17_Mu8Leg2_76X.txt'.format(os.environ['CMSSW_BASE'])
        singleE_path      = '{0}/src/DevTools/Analyzer/data/HWW_HLT_Ele23_WPLoose_76X.txt'.format(os.environ['CMSSW_BASE'])
        doubleELeg1_path  = '{0}/src/DevTools/Analyzer/data/HWW_HLT_Ele17_Ele12Leg1_76X.txt'.format(os.environ['CMSSW_BASE'])
        doubleELeg2_path  = '{0}/src/DevTools/Analyzer/data/HWW_HLT_Ele17_Ele12Leg2_76X.txt'.format(os.environ['CMSSW_BASE'])
        self.hww_singleMu_efficiencies     = self.__parse_hww(singleMu_path,'muons')
        self.hww_doubleMuLeg1_efficiencies = self.__parse_hww(doubleMuLeg1_path,'muons')
        self.hww_doubleMuLeg2_efficiencies = self.__parse_hww(doubleMuLeg2_path,'muons')
        self.hww_singleE_efficiencies      = self.__parse_hww(singleE_path,'electrons')
        self.hww_doubleELeg1_efficiencies  = self.__parse_hww(doubleELeg1_path,'electrons')
        self.hww_doubleELeg2_efficiencies  = self.__parse_hww(doubleELeg2_path,'electrons')
        
        # define supported triggers
        self.singleTriggers = {
            'muons'    : ['IsoMu20_OR_IsoTkMu20', 'Mu45_eta2p1', 'Mu50', 'Mu17_Mu8Leg1', 'Mu17_Mu8Leg2'],
            'electrons': ['Ele23_WPLoose', 'Ele17_Ele12Leg1', 'Ele17_Ele12Leg2'],
            'taus'     : ['PFTau35'],
        }
        self.doubleTriggers = {
            'muons'    : ['Mu17_Mu8'],
            'electrons': ['Ele17_Ele12'],
            'taus'     : ['DoublePFTau35'],
        }

        ################
        ### 80X 2016 ###
        ################
        # tau https://indico.cern.ch/event/544712/contributions/2213574/attachments/1295299/1930984/htt_tau_trigger_17_6_2016.pdf
        self.tau_efficiencies_2016 = {
            'MediumIsoPFTau35_Trk_eta2p1': {
                'val'    : {'m0' : 3.86506E+01, 'sigma' : 5.81155E+00, 'alpha' : 5.82783E+00, 'n' : 3.38903E+00, 'norm' : 9.33449E+00,}, # no iso
            },
            'LooseIsoPFTau20_SingleL1' : {
                'val'    : {'m0' : 2.14111E+01, 'sigma' : 1.05522E+00, 'alpha' : 1.32782E+00, 'n' : 1.50352E+00, 'norm' : 9.96428E-01,}, # no iso
            },
            'LooseIsoPFTau20' : {
                'val'    : {'m0' : 2.13600E+01, 'sigma' : 8.1845E-01, 'alpha' : 6.5401E-01, 'n' : 1.71559E+00, 'norm' : 1.00000E+00,}, # no iso
            },
        }

        # private electron
        self.private_electron_80X = {}
        path = '{0}/src/DevTools/Analyzer/data/trigger_efficiencies_electron_2016.root'.format(os.environ['CMSSW_BASE'])
        self.private_electron_80X_rootfile = ROOT.TFile(path)
        trigMap = {
            'Ele23Ele12DZ'         : 'passingHLTEle23Ele12DZ/probe_sc_pt_probe_sc_eta_PLOT_passingHLTEle23Ele12Leg2_true_&_tag_passingEle23Ele12Leg1_true',
            'Ele23Ele12Leg1'       : 'passingHLTEle23Ele12Leg1/probe_sc_pt_probe_sc_eta_PLOT',
            'Ele23Ele12Leg2'       : 'passingHLTEle23Ele12Leg2/probe_sc_pt_probe_sc_eta_PLOT_tag_passingEle23Ele12Leg1_true',
            'Ele24Tau20LegSingleL1': 'passingHLTEle24Tau20LegSingleL1/probe_sc_pt_probe_sc_eta_PLOT',
            'Ele25Eta2p1'          : 'passingHLTEle25Eta2p1/probe_sc_pt_probe_sc_eta_PLOT',
            'Ele25Tight'           : 'passingHLTEle25Tight/probe_sc_pt_probe_sc_eta_PLOT',
            'Ele35'                : 'passingHLTEle35/probe_sc_pt_probe_sc_eta_PLOT',
        }
        for trig in trigMap:
            self.private_electron_80X[trig] = self.private_electron_80X_rootfile.Get(trigMap[trig])

        # private muon
        self.private_muon_80X = {}
        path = '{0}/src/DevTools/Analyzer/data/trigger_efficiencies_muon_2016.root'.format(os.environ['CMSSW_BASE'])
        self.private_muon_80X_rootfile = ROOT.TFile(path)
        trigMap = {
            'IsoMu22ORIsoTkMu22'  : 'passingIsoMu22ORIsoTkMu22/probe_pt_probe_eta_PLOT',
            'Mu17Mu8Leg1'         : 'passingMu17/probe_pt_probe_eta_PLOT',
            'Mu17Mu8Leg2'         : 'passingMu8/probe_pt_probe_eta_PLOT_tag_passingMu17_true',
            'Mu19Tau20LegSingleL1': 'passingMu19Tau20MLegSingleL1/probe_pt_probe_eta_PLOT',
        }
        for trig in trigMap:
            self.private_muon_80X[trig] = self.private_muon_80X_rootfile.Get(trigMap[trig])


        # define supported triggers
        if self.version=='80X':
            self.singleTriggers = {
                'muons'    : ['IsoMu22ORIsoTkMu22'],
                'electrons': ['Ele25Eta2p1','Ele25Tight','Ele35'],
                'taus'     : [],
            }
            self.doubleTriggers = {
                'muons'    : ['Mu17Mu8','Mu19Tau20SingleL1'],
                'electrons': ['Ele23Ele12','Ele24Tau20SingleL1'],
                'taus'     : ['DoublePFTau35','Ele24Tau20SingleL1','Mu19Tau20SingleL1'],
            }

    def __parse_hww(self,filename,fileType):
        '''Parse text file of trigger efficiencies
           Format:
               muons:     etamin etamax ptmin ptmax eff errup errdown
               electrons: etamin etamax ptmin ptmax eff err
        '''
        scales = []
        with open(filename) as f:
            for line in f.readlines():
                if fileType=='muons':
                    etamin, etamax, ptmin, ptmax, eff, errup, errdown = line.split()
                elif fileType=='electrons':
                    etamin, etamax, ptmin, ptmax, eff, err = line.split()
                    errup = err
                    errdown = err
                else:
                    logging.error('Unrecognized HWW scale fileType: {0}'.format(fileType))
                    return []
                scales += [{
                    'etamin' : float(etamin), 
                    'etamax' : float(etamax), 
                    'ptmin'  : float(ptmin), 
                    'ptmax'  : float(ptmax), 
                    'eff'    : float(eff), 
                    'errup'  : float(errup), 
                    'errdown': float(errdown),
                }]
        return scales

    def __crystalball_fit(self,pt,mode):
        if self.version=='80X' and mode in self.tau_efficiencies_2016: # 2016
            m0    = self.tau_efficiencies_2016[mode]['val']['m0']
            sigma = self.tau_efficiencies_2016[mode]['val']['sigma']
            alpha = self.tau_efficiencies_2016[mode]['val']['alpha']
            n     = self.tau_efficiencies_2016[mode]['val']['n']
            norm  = self.tau_efficiencies_2016[mode]['val']['norm']
        elif self.version=='76X' and mode in self.doubleTau_efficiencies:
            m0    = self.doubleTau_efficiencies[mode]['val']['m0']
            sigma = self.doubleTau_efficiencies[mode]['val']['sigma']
            alpha = self.doubleTau_efficiencies[mode]['val']['alpha']
            n     = self.doubleTau_efficiencies[mode]['val']['n']
            norm  = self.doubleTau_efficiencies[mode]['val']['norm']
        else:
            return 0.
        x = pt
        # recreate the fit
        sqrtPiOver2 = math.sqrt(ROOT.TMath.PiOver2())
        sqrt2 = math.sqrt(2.)
        sig = abs(sigma)
        t = (x-m0)/sig * alpha/abs(alpha)
        absAlpha = abs(alpha/sig)
        a =  ROOT.TMath.Power(n/absAlpha,n) * ROOT.TMath.Exp(-0.5*absAlpha*absAlpha)
        b = absAlpha - n/absAlpha
        arg = absAlpha/sqrt2
        if arg>5.: ApproxErf = 1.
        elif arg<-5.: ApproxErf = -1.
        else: ApproxErf = ROOT.TMath.Erf(arg)
        leftArea = (1.+ApproxErf)*sqrtPiOver2
        rightArea = (a*1./ROOT.TMath.Power(absAlpha-b,n-1))/(n-1)
        area = leftArea+rightArea
        if t<= absAlpha:
            arg = t/sqrt2
            if arg>5.: ApproxErf = 1.
            elif arg<-5.: ApproxErf = -1.
            else: ApproxErf = ROOT.TMath.Erf(arg)
            return norm * (1.+ApproxErf)*sqrtPiOver2/area
        else:
            return norm * (leftArea + a * (1./ROOT.TMath.Power(t-b,n-1) - 1./ROOT.TMath.Power(absAlpha-b,n-1))/(1-n))/area

    def __exit__(self, type, value, traceback):
        self.__finish()

    def __del__(self):
        self.__finish()

    def __finish(self):
        self.singleMu_rootfile.Close()
        self.private_electron_80X_rootfile.Close()
        self.private_muon_80X_rootfile.Close()

    def __triggerWarning(self,triggers):
        logging.warning('Unmatched triggers: {0}'.format(' '.join(triggers)))
        return 0.

    def __getEfficiency(self,rootName,mode,cand):
        pt = cand.pt()
        eta = cand.eta()
        if cand.collName=='muons':
            if self.version=='76X':
                # Muon POG
                # ignore Run2015C, reweight isomu via hlt trigger
                if rootName == 'IsoMu20_OR_IsoTkMu20':
                    if pt>120: pt = 119
                    name0 = 'runD_{0}_HLTv4p2'.format(rootName)
                    name1 = 'runD_{0}_HLTv4p3'.format(rootName)
                    hist0 = self.singleMu_efficiencies[name0][mode]
                    hist1 = self.singleMu_efficiencies[name1][mode]
                    val0 = hist0.GetBinContent(hist0.FindBin(pt,abs(eta)))
                    val1 = hist0.GetBinContent(hist1.FindBin(pt,abs(eta)))
                    return (0.401*val0+1.899*val1)/2.3
                elif rootName in ['Mu50', 'Mu45_eta2p1']:
                    if pt>120: pt = 119
                    hist = self.singleMu_efficiencies['runD_{0}'.format(rootName)][mode]
                    return hist.GetBinContent(hist.FindBin(pt,abs(eta)))
                # HWW
                elif rootName=='Mu17_Mu8Leg1':
                    if pt>200: pt = 199
                    for row in self.hww_doubleMuLeg1_efficiencies:
                       if (eta>=row['etamin'] 
                           and eta<=row['etamax']
                           and pt>=row['ptmin']
                           and pt<=row['ptmax']):
                           return row['eff']
                elif rootName=='Mu17_Mu8Leg2':
                    if pt>200: pt = 199
                    for row in self.hww_doubleMuLeg2_efficiencies:
                       if (eta>=row['etamin'] 
                           and eta<=row['etamax']
                           and pt>=row['ptmin']
                           and pt<=row['ptmax']):
                           return row['eff']
            elif self.version=='80X':
                if rootName in self.private_muon_80X:
                    if pt>1000: pt = 999
                    hist = self.private_muon_80X[rootName]
                    return hist.GetBinContent(hist.FindBin(pt,eta))
        elif cand.collName=='electrons':
            if self.version=='76X':
                # HWW
                if rootName=='Ele23_WPLoose':
                    if pt>100: pt = 99
                    for row in self.hww_singleE_efficiencies:
                       if (eta>=row['etamin'] 
                           and eta<=row['etamax']
                           and pt>=row['ptmin']
                           and pt<=row['ptmax']):
                           return row['eff']
                elif rootName=='Ele17_Ele12Leg1':
                    if pt>100: pt = 99
                    for row in self.hww_doubleELeg1_efficiencies:
                       if (eta>=row['etamin'] 
                           and eta<=row['etamax']
                           and pt>=row['ptmin']
                           and pt<=row['ptmax']):
                           return row['eff']
                elif rootName=='Ele17_Ele12Leg2':
                    if pt>100: pt = 99
                    for row in self.hww_doubleELeg2_efficiencies:
                       if (eta>=row['etamin'] 
                           and eta<=row['etamax']
                           and pt>=row['ptmin']
                           and pt<=row['ptmax']):
                           return row['eff']
            elif self.version=='80X':
                if rootName in self.private_electron_80X:
                    if pt>1000: pt = 999
                    hist = self.private_electron_80X[rootName]
                    return hist.GetBinContent(hist.FindBin(pt,eta))
        elif cand.collName=='taus':
            if self.version=='76X':
                return self.__crystalball_fit(pt,mode)
            elif self.version=='80X':
                return self.__crystalball_fit(pt,rootName)
        return 0.

    def __getLeadEfficiency(self,rootNames,mode,cand):
        if cand.collName=='electrons':
            if 'Ele17_Ele12' in rootNames:
                return self.__getEfficiency('Ele17_Ele12Leg1',mode,cand)
            if 'Ele23Ele12' in rootNames:
                return self.__getEfficiency('Ele23Ele12Leg1',mode,cand)
            if 'Ele24Tau20SingleL1' in rootNames:
                return self.__getEfficiency('Ele24Tau20LegSingleL1',mode,cand)
        elif cand.collName=='muons':
            if 'Mu17_Mu8' in rootNames:
                return self.__getEfficiency('Mu17_Mu8Leg1',mode,cand)
            if 'Mu17Mu8' in rootNames:
                return self.__getEfficiency('Mu17Mu8Leg1',mode,cand)
            if 'Mu19Tau20SingleL1' in rootNames:
                return self.__getEfficiency('Mu19Tau20LegSingleL1',mode,cand)
        elif cand.collName=='taus':
            if 'DoublePFTau35' in rootNames:
                return self.__getEfficiency('MediumIsoPFTau35_Trk_eta2p1',mode,cand)
        return 0.

    def __getTrailEfficiency(self,rootNames,mode,cand):
        if cand.collName=='electrons':
            if 'Ele17_Ele12' in rootNames:
                return self.__getEfficiency('Ele17_Ele12Leg2',mode,cand)
            if 'Ele23Ele12' in rootNames:
                return self.__getEfficiency('Ele23Ele12Leg2',mode,cand)
        elif cand.collName=='muons':
            if 'Mu17_Mu8' in rootNames:
                return self.__getEfficiency('Mu17_Mu8Leg2',mode,cand)
            if 'Mu17Mu8' in rootNames:
                return self.__getEfficiency('Mu17Mu8Leg2',mode,cand)
        elif cand.collName=='taus':
            if 'DoublePFTau35' in rootNames:
                return self.__getEfficiency('MediumIsoPFTau35_Trk_eta2p1',mode,cand)
            if 'Ele24Tau20SingleL1' in rootNames:
                return self.__getEfficiency('LooseIsoPFTau20_SingleL1',mode,cand)
            if 'Mu19Tau20SingleL1' in rootNames:
                return self.__getEfficiency('LooseIsoPFTau20_SingleL1',mode,cand)
        return 0.

    def __getSingleEfficiency(self,rootNames,mode,cand):
        if cand.collName=='electrons':
            if 'Ele23_WPLoose' in rootNames:
                return self.__getEfficiency('Ele23_WPLoose',mode,cand)
            elif 'Ele17_Ele12Leg1' in rootNames:
                return self.__getEfficiency('Ele17_Ele12Leg1',mode,cand)
            elif 'Ele17_Ele12Leg2' in rootNames:
                return self.__getEfficiency('Ele17_Ele12Leg2',mode,cand)
            elif 'Ele25Eta2p1' in rootNames:
                return self.__getEfficiency('Ele25Eta2p1',mode,cand)
            elif 'Ele25Tight' in rootNames:
                return self.__getEfficiency('Ele25Tight',mode,cand)
        elif cand.collName=='muons':
            if 'IsoMu20_OR_IsoTkMu20' in rootNames:
                return self.__getEfficiency('IsoMu20_OR_IsoTkMu20',mode,cand)
            elif 'Mu45_eta2p1' in rootNames:
                return self.__getEfficiency('Mu45_eta2p1',mode,cand)
            elif 'Mu50' in rootNames:
                return self.__getEfficiency('Mu50',mode,cand)
            elif 'Mu17_Mu8Leg1' in rootNames:
                return self.__getEfficiency('Mu17_Mu8Leg1',mode,cand)
            elif 'Mu17_Mu8Leg2' in rootNames:
                return self.__getEfficiency('Mu17_Mu8Leg2',mode,cand)
            elif 'IsoMu22ORIsoTkMu22' in rootNames:
                return self.__getEfficiency('IsoMu22ORIsoTkMu22',mode,cand)
        elif cand.collName=='taus':
            return 0.
        return 0.

    def __hasSingle(self,triggers,triggerType):
        for trigger in triggers:
            if trigger in self.singleTriggers[triggerType]: return True
        return False

    def __hasDouble(self,triggers,triggerType):
        for trigger in triggers:
            if trigger in self.doubleTriggers[triggerType]: return True
        return False


    def __getTriggerEfficiency(self,triggers,cands,mode):
        '''Get an efficiency'''

        #######################
        ### Single triggers ###
        #######################
        if ((len(triggers)==1 and (self.__hasSingle(triggers,'electrons') or self.__hasSingle(triggers,'muons'))) or
            (len(triggers)==2 and (self.__hasSingle(triggers,'electrons') and self.__hasSingle(triggers,'muons')))):
            val = 1-product([1-self.__getSingleEfficiency(triggers,mode,cand) for cand in cands])
        

        #######################
        ### Double triggers ###
        #######################
        elif ((len(triggers)==1 and (self.__hasDouble(triggers,'electrons') or self.__hasDouble(triggers,'muons') or   # ee, mm, tt
                  self.__hasDouble(triggers,'taus'))) or
              (len(triggers)==2 and (self.__hasDouble(triggers,'electrons') and self.__hasDouble(triggers,'muons')))): # ee/mm
            val = 1-(
                # none pass lead
                product([1-self.__getLeadEfficiency(triggers,mode,cand) for cand in cands])
                # one pass lead, none pass trail
                +sum([self.__getLeadEfficiency(triggers,mode,lead)
                      *product([1-self.__getTrailEfficiency(triggers,mode,trail) if trail!=lead else 1. for trail in cands if trail.collName==lead.collName]) for lead in cands])
                # TODO: DZ not included ???
                # one pass lead, one pass trail, fail dz
                # one pass lead, two pass trail, both fail dz
            )

        ################################
        ### single + double triggers ###
        ################################
        elif ((len(triggers)==2 and ((self.__hasSingle(triggers,'electrons') and self.__hasDouble(triggers,'electrons')) or # e/ee
                                     (self.__hasSingle(triggers,'muons') and self.__hasDouble(triggers,'muons')))) or       # m/mm
              (len(triggers)==3 and ((self.__hasSingle(triggers,'electrons') and self.__hasSingle(triggers,'muons')         # e/m/tt
                                     and self.__hasDouble(triggers,'taus')))) or
              (len(triggers)==4 and (self.__hasSingle(triggers,'electrons') and self.__hasDouble(triggers,'electrons') and  # e/m/ee/mm (e/m/ee/em/mm)
                                     self.__hasSingle(triggers,'muons') and self.__hasDouble(triggers,'muons'))) or
              (len(triggers)==5 and (self.__hasSingle(triggers,'electrons') and self.__hasDouble(triggers,'electrons') and  # e/m/ee/mm/tt (e/m/ee/em/mm/tt)
                                     self.__hasSingle(triggers,'muons') and self.__hasDouble(triggers,'muons') and
                                     self.__hasDouble(triggers,'taus')))):
            val = 1-(
                # none pass single
                product([1-self.__getSingleEfficiency(triggers,mode,cand) for cand in cands if cand.collName in ['electrons','muons']]) # only electron/muon single triggers
            )*(
                # none pass lead
                product([1-self.__getLeadEfficiency(triggers,mode,cand) for cand in cands])
                # one pass lead, none pass trail
                +sum([self.__getLeadEfficiency(triggers,mode,lead)
                      *product([1-self.__getTrailEfficiency(triggers,mode,trail) if trail!=lead else 1. for trail in cands if (
                                       (trail.collName=='taus' and lead.collName=='taus') or                                # no cross trigger with taus
                                       (trail.collName in ['electrons','muons'] and lead.collName in ['electrons','muons']) # allow cross triggers with e/m
                                   )
                               ]) for lead in cands])
                # TODO: DZ not included ???
                # one pass lead, one pass trail, fail dz
                # one pass lead, two pass trail, both fail dz
            )

        else:
            val = 1

        return val

    def getMCEfficiency(self,triggers,cands):
        '''Get the efficiency for a set of triggers for a list of candidates in MC'''
        if self.version=='80X': return 1
        return self.__getTriggerEfficiency(triggers,cands,'MC')

    def getDataEfficiency(self,triggers,cands):
        '''Get the efficiency for a set of triggers for a list of candidates in DATA'''
        return self.__getTriggerEfficiency(triggers,cands,'DATA')

    def getRatio(self,triggers,cands):
        '''Get the scale to apply to MC (eff_data/eff_mc)'''
        if self.version=='80X': return 1
        eff_data = self.getDataEfficiency(triggers,cands)
        eff_mc = self.getMCEfficiency(triggers,cands)
        if eff_mc: return eff_data/eff_mc
        return 1.
