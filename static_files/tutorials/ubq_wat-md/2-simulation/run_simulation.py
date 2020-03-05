from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os

# Argument parser

description = """
Runs an MD simulation in explicit solvent.

Loads a .PDB file from from ../../1-model-water/[system_name].pdb.
Stores simulation results in the local directory.

When running on XSEDE, it is not recommended to use this script directly.
Instead, use the submit_simulation.sh script,
which will create a job file that executes this script.
"""

import argparse
parser = argparse.ArgumentParser(\
  description = description)
parser.add_argument('--system_name', \
  default='1ubq', \
  help='The name of the system.')
parser.add_argument('--simulation_time',
  type=int, \
  default=100, \
  help='The amount of time to simulate, in nanoseconds.')
args = parser.parse_args()

# Input files

pdb = app.PDBFile('../../1-model_water/'+args.system_name+'.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# System Configuration

nonbondedMethod = app.PME
nonbondedCutoff = 1.6*unit.nanometers
switchDistance = 1.2*unit.nanometers
constraints = app.HBonds
rigidWater = True
constraintTolerance = 0.0001

# Integration Options

dt = 0.004*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1/unit.picosecond
pressure = 1*unit.atmospheres
barostatInterval = 1000

# Simulation Options

run_number = 0
while os.path.isfile('%s%03d.dcd'%(args.system_name,run_number)):
  run_number += 1
run_name = '%s%03d'%(args.system_name,run_number)
checkpoint_FN = '%s.chk'%(args.system_name)
print('Data will be stored with the prefix '+run_name)

steps = int(args.simulation_time*unit.nanoseconds/dt)
print('Simulation will run for %d steps\n'%steps)

if os.path.isdir('/Users/'):
  # On my laptop use CPU
  platform = mm.Platform.getPlatformByName('CPU')
  properties = {}
else:
  # Otherwise use CUDA
  platform = mm.Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'mixed'}

dcdReporter = app.DCDReporter(run_name+'.dcd', 1000)
dataReporter = app.StateDataReporter(run_name+'.log', 1000, \
  totalSteps=steps, step=True, time=True, speed=True, progress=True, \
  elapsedTime=True, remainingTime=True, \
  potentialEnergy=True, kineticEnergy=True, totalEnergy=True, \
  temperature=True, volume=True, density=True, separator='\t')
checkpointReporter = app.CheckpointReporter(checkpoint_FN, 10000)

# Prepare the Simulation

print('Building system...')
topology = pdb.topology
positions = pdb.positions
system = forcefield.createSystem(topology, \
  nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, \
  switchDistance=switchDistance, \
  hydrogenMass=4*unit.amu, \
  constraints=constraints, rigidWater=rigidWater)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = mm.LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = app.Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)

# If there is a checkpoint, load it

if not os.path.isfile(checkpoint_FN):
  print('Performing energy minimization...')
  simulation.minimizeEnergy()
else:
  print('Loading from checkpoint')
  with open(checkpoint_FN, 'rb') as f:
    simulation.context.loadCheckpoint(f.read())
simulation.currentStep = 0

# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.step(steps)
