import os
import sys
import logging

import ROOT

import operator


class TriggerPrescales(object):
    '''Class to access the trigger prescales for a given trigger.'''

    def __init__(self,version):
        self.version = version
        self.prescales = {}
        self.prescales['76X'] = {
            # HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg2': 2263.552/4.174,
            # HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg1': 2263.552/45.941,
            # HLT_Mu8_TrkIsoVVL_v*
            'Mu17_Mu8Leg2': 2263.552/1.330,
            # HLT_Mu17_TrkIsoVVL_v*
            'Mu17_Mu8Leg1': 2263.552/197.362,
        }
        self.prescales['80X'] = {
            # HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg2': 3993.88/1.851,
            # HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
            'Ele17_Ele12Leg1': 3993.88/16.246,
            # HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*
            'Ele23_Ele12Leg1': 3993.88/1.851,
            # HLT_Mu8_TrkIsoVVL_v*
            'Mu17_Mu8Leg2':    3993.88/4.111,
            # HLT_Mu17_TrkIsoVVL_v*
            'Mu17_Mu8Leg1':    3993.88/84.638,
        }

    def getPrescale(self,trigger):
        if trigger in self.prescales[self.version]:
            return self.prescales[self.version][trigger]
        else:
            return 1.
