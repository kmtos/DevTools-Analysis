import os
import sys
import logging

import ROOT

import operator

from DevTools.Plotter.utilities import getLumi

class TriggerPrescales(object):
    '''Class to access the trigger prescales for a given trigger.'''

    def __init__(self,version):
        self.version = version
        self.prescales = {}
        self.lumi = {
            '76X': getLumi(version='76X'),
            '80X': getLumi(version='80X'),
        }
        self.prescales['76X'] = {
            # HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg2':   4.174,
            # HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg1':  45.941,
            # HLT_Mu8_TrkIsoVVL_v*
            'Mu17_Mu8Leg2'   :   1.330,
            # HLT_Mu17_TrkIsoVVL_v*
            'Mu17_Mu8Leg1'   : 197.362,
        }
        self.prescales['80X'] = {
            # HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg2':   6.162,
            'Ele23_Ele12Leg2':   6.162,
            # HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg1':  30.397,
            # HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*
            'Ele23_Ele12Leg1':  16.430,
            # HLT_Mu8_TrkIsoVVL_v*
            'Mu17_Mu8Leg2'   :   7.832,
            # HLT_Mu17_TrkIsoVVL_v*
            'Mu17_Mu8Leg1'   : 217.553,
        }

    def getPrescale(self,trigger):
        if trigger in self.prescales[self.version]:
            return self.lumi[self.version]/self.prescales[self.version][trigger]
        else:
            return 1.
