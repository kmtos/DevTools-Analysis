# common leptons ids

from Candidates import Electron, Muon, Tau

##################
### HZZ 4l IDs ###
##################
def passHZZLooseElectron(electron):
    if electron.pt()<=7: return False
    if abs(electron.eta())>=2.5: return False
    if abs(electron.dxy())>=0.5: return False
    if abs(electron.dz())>=1.: return False
    return True

def passHZZLooseMuon(muon):
    if muon.pt()<=5: return False
    if abs(muon.eta())>=2.4: return False
    if abs(muon.dxy())>=0.5: return False
    if abs(muon.dz())>=1.: return False
    if not (muon.isGlobalMuon() or
        (muon.isTrackerMuon() and muon.matchedStations()>0)): return False
    if muon.muonBestTrackType()==2: return False
    return True

def passHZZLoose(cand):
    if cand.__class__.__name__=='Electron': return passHZZLooseElectron(cand)
    if cand.__class__.__name__=='Muon':     return passHZZLooseMuon(cand)
    if cand.__class__.__name__=='Tau':      return passHppLooseTau(cand)
    return False

def passHZZTightElectron(electron):
    if not passHZZLooseElectron(electron): return False
    if electron.relPFIsoRhoR03() >= 0.35: return False
    mva = electron.mvaNonTrigValues()
    eta = electron.superClusterEta()
    pt = electron.pt()
    if pt<=10:
        if eta<0.8:
            return mva > -0.265
        elif eta<1.479:
            return mva > -0.556
        else:
            return mva > -0.551
    else:
        if eta<0.8:
            return mva > -0.072
        elif eta<1.479:
            return mva > -0.286
        else:
            return mva > -0.267

def passHZZTightMuon(muon):
    if not passHZZLooseMuon(muon): return False
    if muon.relPFIsoDeltaBetaR03()>=0.35: return False
    return muon.isPFMuon()

def passHZZTight(cand):
    if cand.__class__.__name__=='Electron': return passHZZTightElectron(cand)
    if cand.__class__.__name__=='Muon':     return passHZZTightMuon(cand)
    if cand.__class__.__name__=='Tau':      return passHppMediumTau(cand)
    return False

#########################
### SMP-16-002 WZ ids ###
#########################
def passWZLooseElectron(electron,version='80X'):
    if electron.pt()<=10: return False
    if abs(electron.eta())>=2.5: return False
    #if electron.wwLoose()<0.5: return False
    if electron.cutBasedHLTPreselection()<0.5: return False
    if abs(electron.superClusterEta())<1.479:
        if electron.dxy()>=0.05: return False
        if electron.dz()>=0.10: return False
    else:
        if electron.dxy()>=0.10: return False
        if electron.dz()>=0.20: return False
    return True

def passWZLooseMuon(muon,version='80X'):
    pt = muon.pt()
    if pt<=10: return False
    if abs(muon.eta())>=2.4: return False
    if version=='76X':
        if muon.isMediumMuon()<0.5: return False
    else:
        if muon.isMediumMuonICHEP()<0.5: return False
    if muon.trackIso()/pt>=0.4: return False
    if muon.relPFIsoDeltaBetaR04()>=0.25: return False
    if abs(muon.dz())>=0.1: return False
    pt = muon.pt()
    dxy = muon.dB2D()
    if abs(dxy)>=0.01 and pt<20: return False
    if abs(dxy)>=0.02 and pt>=20: return False
    return True

def passWZLoose(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZLooseElectron(cand,version=version)
    if cand.__class__.__name__=='Muon':     return passWZLooseMuon(cand,version=version)
    return False

def passWZMediumElectron(electron,version='80X'):
    if not passWZLooseElectron(electron,version=version): return False
    if electron.cutBasedMedium()<0.5: return False
    return True

def passWZMediumMuon(muon,version='80X'):
    if not passWZLooseMuon(muon,version=version): return False
    if muon.relPFIsoDeltaBetaR04()>=0.15: return False
    return True

def passWZMedium(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZMediumElectron(cand,version=version)
    if cand.__class__.__name__=='Muon':     return passWZMediumMuon(cand,version=version)
    return False

def passWZTightElectron(electron,version='80X'):
    if not passWZLooseElectron(electron,version=version): return False
    if electron.cutBasedTight()<0.5: return False
    return True

def passWZTightMuon(muon,version='80X'):
    return passWZMediumMuon(muon,version=version)

def passWZTight(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZTightElectron(cand,version=version)
    if cand.__class__.__name__=='Muon':     return passWZTightMuon(cand,version=version)
    return False

###############
### H++ ids ###
###############
def passHppLooseElectron(electron):
    pt = electron.pt()
    if pt<=10: return False
    if electron.cutBasedLoose()<0.5: return False
    eta = electron.eta()
    dxy = electron.dxy()
    dz = electron.dz()
    if abs(eta)<1.479:
        if dxy>=0.05: return False
        if dz>=0.10: return False
    else:
        if dxy>=0.10: return False
        if dz>=0.20: return False
    return True

def passHppLooseMuon(muon,version='80X'):
    pt = muon.pt()
    if pt<=10: return False
    if abs(muon.eta())>=2.4: return False
    if version=='76X':
        if muon.isMediumMuon()<0.5: return False
    else:
        if muon.isMediumMuonICHEP()<0.5: return False
    if muon.trackIso()/pt>=0.4: return False
    if abs(muon.dz())>=0.5: return False
    dxy = muon.dB2D()
    if abs(dxy)>=0.2: return False
    # previous
    #if abs(muon.dz())>=0.1: return False
    #dxy = muon.dB2D()
    #if abs(dxy)>=0.01 and pt<20: return False
    #if abs(dxy)>=0.02 and pt>=20: return False
    #if muon.relPFIsoDeltaBetaR04()>=0.25: return False
    if muon.relPFIsoDeltaBetaR04()>=0.4: return False
    return True

def passHppLooseTau(tau,new=False):
    pt = tau.pt()
    if pt<=20: return False
    if new:
        if tau.decayModeFindingNewDMs()<0.5: return False
    else:
        if tau.decayModeFinding()<0.5: return False
    if tau.againstMuonLoose3()<0.5: return False
    if tau.againstElectronVLooseMVA6()<0.5: return False
    if abs(tau.dxy())>=0.2: return False
    if abs(tau.dz())>=0.5: return False
    # remove iso
    #if tau.byVLooseIsolationMVArun2v1DBoldDMwLT()<0.5: return False
    #if tau.byIsolationMVArun2v1DBoldDMwLTraw()<-0.8: return False # custom vvloose isolation
    if new:
        if tau.byIsolationMVArun2v1DBnewDMwLTraw()<-0.2: return False # custom vvloose isolation
    else:
        if tau.byIsolationMVArun2v1DBoldDMwLTraw()<-0.2: return False # custom vvloose isolation
    return True

def passHppLoose(cand,version='80X',new=False):
    if cand.__class__.__name__=='Electron': return passHppLooseElectron(cand)
    if cand.__class__.__name__=='Muon':     return passHppLooseMuon(cand,version=version) #if version=='76X' else passHppLooseMuonICHEP(cand)
    if cand.__class__.__name__=='Tau':      return passHppLooseTau(cand,new=new)
    return False

def passHppMediumElectron(electron):
    if not passHppLooseElectron(electron): return False
    return electron.cutBasedMedium()>0.5

def passHppMediumMuon(muon,version='80X'):
    if not passHppLooseMuon(muon,version=version): return False
    if muon.relPFIsoDeltaBetaR04()>=0.15: return False
    return True

def passHppMediumTau(tau):
    if tau.decayModeFinding()<0.5: return False
    if tau.againstMuonLoose3()<0.5: return False
    if tau.againstElectronVLooseMVA6()<0.5: return False
    if tau.byLooseIsolationMVArun2v1DBoldDMwLT()<0.5: return False
    return True

def passHppMedium(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passHppMediumElectron(cand)
    if cand.__class__.__name__=='Muon':     return passHppMediumMuon(cand,version=version) #if version=='76X' else passHppMediumMuonICHEP(cand)
    if cand.__class__.__name__=='Tau':      return passHppMediumTau(cand)
    return False

def passHppTightElectron(electron):
    if not passHppLooseElectron(electron): return False
    return electron.cutBasedTight()>0.5

def passHppTightMuon(muon,version='80X'):
    if not passHppMediumMuon(muon,version=version): return False
    if muon.isTightMuon()<0.5: return False
    return True

def passHppTightTau(tau):
    pt = tau.pt()
    if pt<=20: return False
    if tau.decayModeFinding()<0.5: return False
    if tau.againstMuonLoose3()<0.5: return False
    if tau.againstElectronVLooseMVA6()<0.5: return False
    if tau.byVTightIsolationMVArun2v1DBoldDMwLT()<0.5: return False
    return True

def passHppTight(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passHppTightElectron(cand)
    if cand.__class__.__name__=='Muon':     return passHppTightMuon(cand,version=version) #if version=='76X' else passHppMediumMuonICHEP(cand)
    if cand.__class__.__name__=='Tau':      return passHppTightTau(cand)
    return False

###############
### WZ 2017 ###
###############
def passWZ2017LooseElectron(electron):
    if electron.pt()<=10: return False
    if abs(electron.eta())>=2.5: return False
    if electron.cutBasedVeto()<0.5: return False
    dxy = electron.dxy()
    dz = electron.dz()
    if abs(electron.superClusterEta())<1.479:
        if abs(dxy)>=0.05: return False
        if abs(dz)>=0.10: return False
    else:
        if abs(dxy)>=0.10: return False
        if abs(dz)>=0.20: return False
    return True

def passWZ2017LooseMuon(muon):
    pt = muon.pt()
    if pt<=10: return False
    if abs(muon.eta())>=2.4: return False
    if muon.isMediumMuonICHEP()<0.5: return False
    if muon.trackIso()/pt>=0.4: return False
    if abs(muon.dz())>=0.1: return False
    pt = muon.pt()
    #dxy = muon.dB2D()
    dxy = muon.dxy()
    if abs(dxy)>=0.01 and pt<20: return False
    if abs(dxy)>=0.02 and pt>=20: return False
    if muon.relPFIsoDeltaBetaR04()>=0.25: return False
    return True

def passWZ2017Loose(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZ2017LooseElectron(cand)
    if cand.__class__.__name__=='Muon':     return passWZ2017LooseMuon(cand)
    return False


def passWZ2017MediumElectron(electron):
    if not passWZ2017LooseElectron(electron): return False
    if electron.cutBasedMedium()<0.5: return False
    return True

def passWZ2017MediumMuon(muon):
    if not passWZ2017LooseMuon(muon): return False
    if muon.relPFIsoDeltaBetaR04()>=0.15: return False
    return True

def passWZ2017Medium(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZ2017MediumElectron(cand)
    if cand.__class__.__name__=='Muon':     return passWZ2017MediumMuon(cand)
    return False


def passWZ2017TightElectron(electron):
    if not passWZ2017MediumElectron(electron): return False
    if electron.cutBasedTight()<0.5: return False
    return True

def passWZ2017TightMuon(muon):
    return passWZ2017MediumMuon(muon)

def passWZ2017Tight(cand,version='80X'):
    if cand.__class__.__name__=='Electron': return passWZ2017TightElectron(cand)
    if cand.__class__.__name__=='Muon':     return passWZ2017TightMuon(cand)
    return False


