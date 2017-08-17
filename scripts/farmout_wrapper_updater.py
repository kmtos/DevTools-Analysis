#!/usr/bin/env python
import argparse
import sys
import subprocess
import os

# import run script
from DevTools.Analyzer.runUpdater import runUpdater

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit analyzers')

    parser.add_argument('analysis', type=str, default='', help='Analysis to submit')

    return parser.parse_args(argv)



def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # add on the input and output
    argv = ['--inputFileList',os.environ['INPUT'],'--outputFile',os.environ['OUTPUT']]

    # run the analyzer
    status = runUpdater(args.analysis,argv)

    return status

if __name__ == "__main__":
    func = main()
    sys.exit(func)
