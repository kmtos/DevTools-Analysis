import ROOT

from utilities import deltaR, deltaPhi

##############################
### Basic candidate access ###
##############################
class Candidate(object):
    '''
    Encapsulate access to an object in a TTree.
    '''
    def __init__(self,tree,entry=-1,collName=''):
        self.tree = tree
        self.collName = collName
        self.entry = entry

    def __getattr__(self,name):
        return lambda: self.get(name) # returns the attribute as a function

    def get(self,var):
        '''Default variable access from tree.'''
        varName = '{0}_{1}'.format(self.collName,var)
        return getattr(self.tree,varName)[self.entry]

    def p4(self):
        return ROOT.TLorentzVector(self.pt(), self.eta(), self.phi(), self.energy())

############
### Muon ###
############
class Muon(Candidate):
    '''
    Muon object access.
    '''
    def __init__(self,tree,entry,collName='muons',shift=None):
        super(Muon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def pt(self):
        var = 'rochesterPt'
        if self.shift=='MuonEnUp': var = 'pt_muonEnUp'
        if self.shift=='MuonEnDown': var = 'pt_muonEnDown'
        return self.get(var)

    def eta(self):
        var = 'rochesterEta'
        if self.shift=='MuonEnUp': var = 'eta_muonEnUp'
        if self.shift=='MuonEnDown': var = 'eta_muonEnDown'
        return self.get(var)

    def phi(self):
        var = 'rochesterPhi'
        if self.shift=='MuonEnUp': var = 'phi_muonEnUp'
        if self.shift=='MuonEnDown': var = 'phi_muonEnDown'
        return self.get(var)

    def energy(self):
        var = 'rochesterEnergy'
        if self.shift=='MuonEnUp': var = 'energy_muonEnUp'
        if self.shift=='MuonEnDown': var = 'energy_muonEnDown'
        return self.get(var)

################
### Electron ###
################
class Electron(Candidate):
    '''
    Electron object access.
    '''
    def __init__(self,tree,entry,collName='electrons',shift=None):
        super(Electron, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def pt(self):
        var = 'pt'
        if self.shift=='ElectronEnUp': var = 'pt_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'pt_electronEnDown'
        return self.get(var)

    def eta(self):
        var = 'eta'
        if self.shift=='ElectronEnUp': var = 'eta_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'eta_electronEnDown'
        return self.get(var)

    def phi(self):
        var = 'phi'
        if self.shift=='ElectronEnUp': var = 'phi_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'phi_electronEnDown'
        return self.get(var)

    def energy(self):
        var = 'energy'
        if self.shift=='ElectronEnUp': var = 'energy_electronEnUp'
        if self.shift=='ElectronEnDown': var = 'energy_electronEnDown'
        return self.get(var)

###########
### Tau ###
###########
class Tau(Candidate):
    '''
    Tau object access.
    '''
    def __init__(self,tree,entry,collName='taus',shift=None):
        super(Tau, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def pt(self):
        var = 'pt'
        if self.shift=='TauEnUp': var = 'pt_tauEnUp'
        if self.shift=='TauEnDown': var = 'pt_tauEnDown'
        return self.get(var)

    def eta(self):
        var = 'eta'
        if self.shift=='TauEnUp': var = 'eta_tauEnUp'
        if self.shift=='TauEnDown': var = 'eta_tauEnDown'
        return self.get(var)

    def phi(self):
        var = 'phi'
        if self.shift=='TauEnUp': var = 'phi_tauEnUp'
        if self.shift=='TauEnDown': var = 'phi_tauEnDown'
        return self.get(var)

    def energy(self):
        var = 'energy'
        if self.shift=='TauEnUp': var = 'energy_tauEnUp'
        if self.shift=='TauEnDown': var = 'energy_tauEnDown'
        return self.get(var)

###########
### Jet ###
###########
class Jet(Candidate):
    '''
    Jet object access.
    '''
    def __init__(self,tree,entry,collName='jets',shift=None):
        super(Jet, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def pt(self):
        var = 'pt'
        if self.shift=='JetEnUp': var = 'pt_jetEnUp'
        if self.shift=='JetEnDown': var = 'pt_jetEnDown'
        return self.get(var)

    def eta(self):
        var = 'eta'
        if self.shift=='JetEnUp': var = 'eta_jetEnUp'
        if self.shift=='JetEnDown': var = 'eta_jetEnDown'
        return self.get(var)

    def phi(self):
        var = 'phi'
        if self.shift=='JetEnUp': var = 'phi_jetEnUp'
        if self.shift=='JetEnDown': var = 'phi_jetEnDown'
        return self.get(var)

    def energy(self):
        var = 'energy'
        if self.shift=='JetEnUp': var = 'energy_jetEnUp'
        if self.shift=='JetEnDown': var = 'energy_jetEnDown'
        return self.get(var)

##############
### Photon ###
##############
class Photon(Candidate):
    '''
    Photon object access.
    '''
    def __init__(self,tree,entry,collName='photons',shift=None):
        super(Photon, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

###########
### MET ###
###########
class Met(Candidate):
    '''
    Met object access.
    '''
    def __init__(self,tree,entry=0,collName='pfmet',shift=None):
        super(Met, self).__init__(tree,entry=entry,collName=collName)
        self.shift = shift

    def et(self):
        var = 'et'
        if self.shift=='ElectronEnUp':      var = 'et_electronEnUp'
        if self.shift=='ElectronEnDown':    var = 'et_electronEnDown'
        if self.shift=='MuonEnUp':          var = 'et_muonEnUp'
        if self.shift=='MuonEnDown':        var = 'et_muonEnDown'
        if self.shift=='TauEnUp':           var = 'et_tauEnUp'
        if self.shift=='TauEnDown':         var = 'et_tauEnDown'
        if self.shift=='PhotonEnUp':        var = 'et_photonEnUp'
        if self.shift=='PhotonEnDown':      var = 'et_photonEnDown'
        if self.shift=='JetEnUp':           var = 'et_jetEnUp'
        if self.shift=='JetEnDown':         var = 'et_jetEnDown'
        if self.shift=='JetResUp':          var = 'et_jetResUp'
        if self.shift=='JetResDown':        var = 'et_jetResDown'
        if self.shift=='UnclusteredEnUp':   var = 'et_unclusteredEnUp'
        if self.shift=='UnclusteredEnDown': var = 'et_unclusteredEnDown'
        return self.get(var)

    def phi(self):
        var = 'phi'
        if self.shift=='ElectronEnUp':      var = 'phi_electronEnUp'
        if self.shift=='ElectronEnDown':    var = 'phi_electronEnDown'
        if self.shift=='MuonEnUp':          var = 'phi_muonEnUp'
        if self.shift=='MuonEnDown':        var = 'phi_muonEnDown'
        if self.shift=='TauEnUp':           var = 'phi_tauEnUp'
        if self.shift=='TauEnDown':         var = 'phi_tauEnDown'
        if self.shift=='PhotonEnUp':        var = 'phi_photonEnUp'
        if self.shift=='PhotonEnDown':      var = 'phi_photonEnDown'
        if self.shift=='JetEnUp':           var = 'phi_jetEnUp'
        if self.shift=='JetEnDown':         var = 'phi_jetEnDown'
        if self.shift=='JetResUp':          var = 'phi_jetResUp'
        if self.shift=='JetResDown':        var = 'phi_jetResDown'
        if self.shift=='UnclusteredEnUp':   var = 'phi_unclusteredEnUp'
        if self.shift=='UnclusteredEnDown': var = 'phi_unclusteredEnDown'
        return self.get(var)

    def p4(self):
        metP4 = ROOT.TLorentzVector()
        metP4.SetPtEtaPhiM(self.et(),0.,self.phi(),0)
        return metP4

############################
### Composite candidates ###
############################
class CompositeCandidate(object):
    '''
    Primary object for access to composite variables.
    '''
    def __init__(self,*objects):
        self.objects = objects

    def __getattr__(self,name):
        try:
            return self.get(name)
        except:
            return self.get(name.capitalize()) # catch things like 'pt' instead of 'Pt'

    def get(self,var):
        '''Default variable access from TLorentzVector'''
        vec = self.p4()
        return getattr(vec,var)

    def p4(self):
        return sum([obj.p4() for obj in self.objects])

###################
### Dicandidate ###
###################
class DiCandidate(CompositeCandidate):
    '''
    Dicandidate variable access.
    '''
    def __init__(self,obj0,obj1):
        super(DiCandidate, self).__init__(obj0,obj1)

    def deltaR(self):
        return deltaR(self.objects[0].eta(),
                      self.objects[0].phi(),
                      self.objects[1].eta(),
                      self.objects[1].phi())

    def deltaPhi(self):
        return deltaPhi(self.objects[0].phi(),
                       self.objects[1].phi())

    def deltaEta(self):
        return abs(self.objects[0].eta()-self.objects[1].eta())

##########################
### Candidate plus met ###
##########################
class MetCompositeCandidate(CompositeCandidate):
    '''
    Met + candidate variable specialization.
    '''
    def __init__(self,met,*objects):
        super(MetCompositeCandidate, self).__init__(met,*objects)
        self.met = met
        self.cands = objects

    def metP4(self):
        return self.met.p4()

    def candP4(self):
        return sum([cand.p4() for cand in self.cands])

    def mt(self):
        metP4 = self.metP4()
        candP4 = self.candP4()
        return math.sqrt(abs((candP4.Et()+metP4.Et())**2 - ((candP4+metP4).Pt())**2))

    def Mt(self):
        return self.mt()

    def deltaPhi(self):
        candP4 = self.candP4()
        return deltaPhi(self.met.phi(),candP4.Phi())
