#!/usr/bin/python

import argparse
parser = argparse.ArgumentParser(description = \
  'Create a script or set of scripts to run AutoDock Vina on XSEDE Bridges')
parser.add_argument('--receptor', \
  default='hsg1.pdbqt',
  help='Receptor pdbqt file (Vina argument)')
parser.add_argument('--config', \
  default='config.txt', \
  help='Vina configuration file (Vina argument)')
parser.add_argument('--cpu', \
  default=None,
  type=int,
  help='Number of CPUs to use (Vina argument)')
parser.add_argument('--library', \
  default='library', \
  help='Directory containing library of pdbqt files')
parser.add_argument('--output', \
  default='poses',
  help='Directory containing library of docked poses')
parser.add_argument('--nscripts', \
  type=int, default=1,
  help='Number of scripts, which can be equal to the number of cores')
parser.add_argument('-f', default='Dummy argument')
args = parser.parse_args()

import os
import glob

if args.cpu is None:
  cpu_string = ''
else:
  cpu_string = ' --cpu %d'%args.cpu

lib_FNs = glob.glob(os.path.join(args.library,'*.pdbqt'))
for n in range(args.nscripts):
  sh_F = open('script%d.sh'%n,'w')
  sh_F.write('#!/bin/bash\n')
  sh_F.write('module load autodock/4.2.6\n')
  for lib_FN in lib_FNs[n::args.nscripts]:
    out_FN = os.path.join(args.output, \
      os.path.basename(args.receptor)[:-6] + '_' + \
      os.path.basename(lib_FN)[:-6] + '.pdbqt')
    if not os.path.isfile(out_FN):
      sh_F.write('vina --config ' + args.config + \
        ' --receptor ' + args.receptor + \
        ' --ligand ' + lib_FN + \
        ' --out ' + out_FN + \
        cpu_string + \
        '\n')
