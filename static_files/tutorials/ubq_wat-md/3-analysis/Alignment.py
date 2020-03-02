# Argument parser

import argparse
parser = argparse.ArgumentParser(\
  description = 'Concatenates and aligns a series of simulations in explicit solvent.' + \
    'Loads a .PDB file from from ../../1-model-water/[system_name].pdb.' + \
    'Loads .DCD files from from ../../2-simulation/[system_name][simulation_number].dcd.' + \
    'Stores aligned simulation results in the local directory.')
parser.add_argument('--system_name', \
  default='1ubq', \
  help='The name of the system.')
args = parser.parse_args()

import os

pdb_FN = '../1-model_water/%s.pdb'%args.system_name
aligned_sim_FN = '0/%s_aligned.dcd'%args.system_name
aligned_sim_protein_FN = '0/%s_aligned_protein.dcd'%args.system_name

# Creates list of trajectory files
sim_FNs = []
run_number = 0
sim_FN = '../2-simulation/0/%s%03d.dcd'%(args.system_name,run_number)
while os.path.isfile(sim_FN):
  sim_FNs.append(sim_FN)
  run_number += 1
  sim_FN = '../2-simulation/0/%s%03d.dcd'%(args.system_name,run_number)
print('The following trajectory files will be loaded:', sim_FNs)

# Loads the trajectories
import MDAnalysis as mda
ref = mda.Universe(pdb_FN)
sim = mda.Universe(pdb_FN, sim_FNs)

# Write protein atoms to a new file
ref.select_atoms("protein").write('protein.pdb', frames='all')

# Creates the output directory
import os
if not os.path.isdir('0'):
  os.mkdir('0')

# Align all frames to the reference,
#   minimizing the alpha carbon RMSD
#   and storing the trajectory to a dcd file
from MDAnalysis.analysis import align
alignment = align.AlignTraj(sim, ref, \
  select="protein and name CA", \
  filename=aligned_sim_FN).run()
del sim

# Write protein atoms to a new file
aligned_sim = mda.Universe(pdb_FN, aligned_sim_FN)
aligned_sim.select_atoms("protein").write(aligned_sim_protein_FN, frames='all')
