#!/usr/bin/env python
import os
import sys
import glob
import argparse
from prettytable import PrettyTable


from DevTools.Analyzer.utilities import getNtupleDirectory
from DevTools.Plotter.xsec import getXsec
from DevTools.Utilities.utilities import getCMSSWVersion

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Dump samples available')

    parser.add_argument('--version',type=str,default=getCMSSWVersion(),help='Samples to dump')

    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    table = PrettyTable(['Sample','xsec [pb]'])
    table.align = 'r'
    table.align['Sample'] = 'l'

    for sample in sorted(glob.glob(os.path.join(getNtupleDirectory(version=args.version),'*'))):
        name = os.path.basename(sample)
        table.add_row([name,getXsec(name)])

    print table.get_string()


if __name__ == "__main__":
    status = main()
    sys.exit(status)
