# common leptons ids

from Candidates import Electron, Muon, Tau

#########################
### SMP-16-002 WZ ids ###
#########################
def passWZLooseElectron(electron):
    if electron.pt()<=10: return False
    if abs(electron.eta())>=2.5: return False
    if electron.wwLoose()<0.5: return False
    return True

def passWZLooseMuon(muon):
    pt = muon.pt()
    if pt<=10: return False
    if abs(muon.eta())>=2.4: return False
    if muon.isMediumMuon()<0.5: return False
    if muon.trackIso()/pt>=0.4: return False
    if muon.relPFIsoDeltaBetaR04()>=0.4: return False
    return True

def passWZLoose(cand):
    if isinstance(cand,Electron): return passWZLooseElectron(cand)
    if isinstance(cand,Muon):     return passWZLooseMuon(cand)
    return False

def passWZMediumElectron(electron):
    if not passWZLooseElectron(electron): return False
    if electron.cutBasedMedium()<0.5: return False
    return True

def passWZMediumMuon(muon):
    if not passWZLooseMuon(muon): return False
    if abs(muon.dz())>=0.1: return False
    pt = muon.pt()
    dxy = muon.dB2D()
    if abs(dxy)>=0.01 and pt<20: return False
    if abs(dxy)>=0.02 and pt>=20: return False
    if muon.relPFIsoDeltaBetaR04()>=0.15: return False
    return True

def passWZMedium(cand):
    if isinstance(cand,Electron): return passWZMediumElectron(cand)
    if isinstance(cand,Muon):     return passWZMediumMuon(cand)
    return False

def passWZTightElectron(electron):
    if not passWZLooseElectron(electron): return False
    if electron.cutBasedTight()<0.5: return False
    return True

def passWZTightMuon(muon):
    return passWZMediumMuon(muon)

def passWZTight(cand):
    if isinstance(cand,Electron): return passWZTightElectron(cand)
    if isinstance(cand,Muon):     return passWZTightMuon(cand)
    return False

###############
### H++ ids ###
###############
def passHppLooseElectron(electron):
    return electron.cutBasedVeto()>0.5

def passHppLooseMuon(muon):
    return passWZLooseMuon(muon)

def passHppLooseTau(tau):
        if tau.decayModeFinding()<0.5: return False
        if tau.againstMuonLoose3()<0.5: return False
        if tau.againstElectronVLooseMVA6()<0.5: return False
        # remove iso
        # if tau.byLooseIsolationMVArun2v1DBoldDMwLT()<0.5: return False
        return True

def passHppLoose(cand):
    if isinstance(cand,Electron): return passHppLooseElectron(cand)
    if isinstance(cand,Muon):     return passHppLooseMuon(cand)
    if isinstance(cand,Tau):      return passHppLooseTau(cand)
    return False

def passHppMediumElectron(electron):
    if not passHppLooseElectron(electron): return False
    return electron.cutBasedMedium()>0.5

def passHppMediumMuon(muon):
    if not passHppLooseMuon(muon): return False
    return passWZMediumMuon(muon)

def passHppMediumTau(tau):
    if tau.decayModeFinding()<0.5: return False
    if tau.againstMuonLoose3()<0.5: return False
    if tau.againstElectronVLooseMVA6()<0.5: return False
    if tau.byLooseIsolationMVArun2v1DBoldDMwLT()<0.5: return False
    return True

def passHppMedium(cand):
    if isinstance(cand,Electron): return passHppMediumElectron(cand)
    if isinstance(cand,Muon):     return passHppMediumMuon(cand)
    if isinstance(cand,Tau):      return passHppMediumTau(cand)
    return False

def passHppTightElectron(electron):
    if not passHppLooseElectron(electron): return False
    return electron.cutBasedTight()>0.5

def passHppTightMuon(muon):
    if not passHppLooseMuon(muon): return False
    return passWZMediumMuon(muon)

def passHppTightTau(tau):
    if tau.decayModeFinding()<0.5: return False
    if tau.againstMuonLoose3()<0.5: return False
    if tau.againstElectronVLooseMVA6()<0.5: return False
    if tau.byVTightIsolationMVArun2v1DBoldDMwLT()<0.5: return False
    return True

def passHppTight(cand):
    if isinstance(cand,Electron): return passHppTightElectron(cand)
    if isinstance(cand,Muon):     return passHppTightMuon(cand)
    if isinstance(cand,Tau):      return passHppTightTau(cand)
    return False
