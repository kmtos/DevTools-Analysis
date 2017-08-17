import logging
import os
import sys
import argparse

from array import array

from DevTools.Analyzer.TreeUpdaterBase import TreeUpdaterBase
from DevTools.Plotter.utilities import getTestFiles
from DevTools.Analyzer.Candidates import *

import ROOT

logger = logging.getLogger("EGTreeUpdater")
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

class EGTreeUpdater(TreeUpdaterBase):

    def __init__(self,**kwargs):
        inputTreeName = kwargs.pop('inputTreeName','EGTree')
        outputFileName = kwargs.pop('outputFileName','egTreeUpdated.root')
        outputTreeName = kwargs.pop('outputTreeName','EGTree')
        super(EGTreeUpdater, self).__init__(outputFileName=outputFileName,outputTreeName=outputTreeName,inputTreeName=inputTreeName,**kwargs)

        # setup rows to update
        self.phos = ['g']
        self.scales = ['preselection','mva','mvaPreselection']
        self.scaleNames = {
            'preselection'   : 'Preselection',
            'mva'            : 'MVA0p0',
            'mvaPreselection': 'MVA0p0Pre',
        }

        self.allVars = {}
        self.allBranches = {}
        for pho in self.phos:
            for scale in self.scales:
                for shift in ['','Up','Down']:
                    name = '{0}_{1}Scale{2}'.format(pho,scale,shift)
                    self.allVars[name] = array('f',[0])
                    if hasattr(self.tree,name): # update branch value
                        self.tree.SetBranchAddress(name,self.allVars[name])
                    else: # new branch
                        self.tree.Branch(name, self.allVars[name], "{0}/F".format(name))

    def perRowAction(self):
        # create candidates and update values
        for pho in self.phos:
            g = Photon(self.oldtree,collName=pho)
            for scale in self.scales:
                s = self.leptonScales.getScale(self.scaleNames[scale],g,doError=True)
                for i,shift in enumerate(['','Up','Down']):
                    name = '{0}_{1}Scale{2}'.format(pho,scale,shift)
                    self.allVars[name][0] = s[i]

        # fill tree
        self.tree.Fill()




def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Run analyzer')

    parser.add_argument('--inputFiles', type=str, nargs='*', default=getTestFiles('EG','dy',version='80X'), help='Input files')
    parser.add_argument('--inputFileList', type=str, default='', help='Input file list')
    parser.add_argument('--outputFile', type=str, default='egTreeUpdated.root', help='Output file')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    egAnalysis = EGTreeUpdater(
        outputFileName=args.outputFile,
        outputTreeName='EGTree',
        inputFileNames=args.inputFileList if args.inputFileList else args.inputFiles,
        inputTreeName='EGTree',
    )

    try:
       egAnalysis.update()
    except KeyboardInterrupt:
       egAnalysis.finish()

    return 0

if __name__ == "__main__":
    status = main()
    sys.exit(status)
