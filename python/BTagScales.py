import os
import sys
import logging
import csv

import ROOT

from DevTools.Utilities.utilities import *

class BTagScales(object):
    '''Class to access btag scale factors'''

    def __init__(self,version):
        self.version = version

        ROOT.gSystem.Load('libCondFormatsBTauObjects')
        ROOT.gSystem.Load('libCondToolsBTau')

        self.v_sys = getattr(ROOT, 'vector<string>')()
        self.v_sys.push_back('up')
        self.v_sys.push_back('down')

        # B tag POG
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
        # https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration

        # 80X
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
        if self.version == '80X':
            # CSVv2
            path = '{0}/src/DevTools/Analyzer/data/CSVv2_Moriond17_B_H.csv'.format(os.environ['CMSSW_BASE'])
            self.sf_reader_csv = {}
            for wp in [0,1,2]:
                self.sf_reader_csv[wp] = self._get_reader_csv(path,wp)


    def _get_reader_csv(self,path,wp):
        calib = ROOT.BTagCalibration('csvv2',path)
        reader = ROOT.BTagCalibrationReader(wp,"central",self.v_sys)
        reader.load(calib,0,'comb')
        reader.load(calib,1,'comb')
        reader.load(calib,2,'incl')
        return reader

    def get_sf_csv(self,wp,pt,eta,flavor=0,shift='central'):
        return self.sf_reader_csv[wp].eval_auto_bounds(shift,flavor,pt,eta)


if __name__ == "__main__":
    scales = BTagScales(getCMSSWVersion())
    print scales.get_sf_csv(1,1.2,31)
    print scales.get_sf_csv(1,1.2,31,shift='up')
    print scales.get_sf_csv(1,1.2,31,shift='down')
    print scales.get_sf_csv(1,1.2,31,flavor=1)
    print scales.get_sf_csv(1,1.2,31,flavor=1,shift='up')
    print scales.get_sf_csv(1,1.2,31,flavor=1,shift='down')
    print scales.get_sf_csv(1,1.2,31,flavor=2)
    print scales.get_sf_csv(1,1.2,31,flavor=2,shift='up')
    print scales.get_sf_csv(1,1.2,31,flavor=2,shift='down')
    sys.exit(0)
