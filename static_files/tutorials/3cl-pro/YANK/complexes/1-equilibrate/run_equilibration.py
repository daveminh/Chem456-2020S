from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from parmed.openmm.reporters import RestartReporter
from sys import stdout
import os

# Argument parser

description = """
Runs an MD simulation in explicit solvent.

Loads AMBER files [system_name].prmtop and [system_name].inpcrd
Stores simulation results in the local directory.

When running on XSEDE, it is not recommended to use this script directly.
Instead, use the submit_simulation.sh script,
which will create a job file that executes this script.
"""

import argparse
parser = argparse.ArgumentParser(\
  description = description)
parser.add_argument('--system_name', \
  default='../0-build/ZINC000001542916_solv', \
  help='The base path of the system.')
parser.add_argument('--simulation_time',
  type=int, \
  default=5, \
  help='The amount of time to simulate, in nanoseconds.')
args = parser.parse_args()

if not (os.path.isfile(args.system_name+'.prmtop') and \
    os.path.isfile(args.system_name+'.inpcrd')):
  raise Exception('Input files not found')

# System Configuration

nonbondedMethod = app.PME
nonbondedCutoff = 1.0*unit.nanometers
switchDistance = 0.9*unit.nanometers
constraints = app.HBonds
rigidWater = True
constraintTolerance = 0.0001

# Integration Options

dt = 0.002*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1/unit.picosecond
pressure = 1*unit.atmospheres
barostatInterval = 1000

# Simulation Options

run_number = 0
while os.path.isfile('%s%03d.dcd'%(os.path.basename(args.system_name),run_number)):
  run_number += 1
run_name = '%s%03d'%(os.path.basename(args.system_name),run_number)
checkpoint_FN = '%s.chk'%(os.path.basename(args.system_name))
rst7_FN = '%s.rst7'%(os.path.basename(args.system_name))
print('Data will be stored with the prefix '+run_name)

steps = int(args.simulation_time*unit.nanoseconds/dt)
print('Simulation will run for %d steps\n'%steps)

if os.path.isdir('/Users/'):
  # On my laptop use OpenCL
  platform = mm.Platform.getPlatformByName('OpenCL')
  properties = {'Precision': 'mixed'}
else:
  # Otherwise use CUDA
  platform = mm.Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'mixed'}

# Report
#   dcd every 10 ps
#   state every 10 ps
#   checkpoint every 100 ps
dcdReporter = app.DCDReporter(run_name+'.dcd', 5000)
dataReporter = app.StateDataReporter(run_name+'.log', 5000, \
  totalSteps=steps, step=True, time=True, speed=True, progress=True, \
  elapsedTime=True, remainingTime=True, \
  potentialEnergy=True, kineticEnergy=True, totalEnergy=True, \
  temperature=True, volume=True, density=True, separator='\t')
checkpointReporter = app.CheckpointReporter(checkpoint_FN, 50000)
restartReporter = RestartReporter(rst7_FN, 50000)

# Prepare the Simulation

print('Setting up system...')
prmtop = app.AmberPrmtopFile(args.system_name+'.prmtop')
inpcrd = app.AmberInpcrdFile(args.system_name+'.inpcrd')

system = prmtop.createSystem(\
  nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, \
  switchDistance=switchDistance, \
  hydrogenMass=4*unit.amu, \
  constraints=constraints, rigidWater=rigidWater)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = mm.LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

# If there is a checkpoint, load it

if not os.path.isfile(checkpoint_FN):
  print('Performing energy minimization...')
  state = simulation.context.getState(getEnergy=True, getPositions=True)
  Eo = state.getPotentialEnergy()
  simulation.minimizeEnergy()
  state = simulation.context.getState(getEnergy=True, getPositions=True)
  En = state.getPotentialEnergy()
  print(f'Minimized energy from {Eo} to {En}')
  positions = state.getPositions()

  print('Ramping temperature')
  for lT in ([25] + [t for t in range(50,300,50)]):
    system_lT = prmtop.createSystem(\
      nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, \
      switchDistance=switchDistance, \
      hydrogenMass=4*unit.amu, \
      constraints=constraints, rigidWater=rigidWater)
    system_lT.addForce(mm.MonteCarloBarostat(pressure, lT * unit.kelvin, barostatInterval))
    integrator_lT = mm.LangevinIntegrator(lT * unit.kelvin, friction, dt)
    integrator_lT.setConstraintTolerance(constraintTolerance)
    simulation_lT = app.Simulation(prmtop.topology, system_lT, integrator_lT, platform, properties)
    simulation_lT.context.setPositions(positions)
    simulation_lT.step(5000)

    state = simulation_lT.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
    print('After equilibration at %d K, potential energy is'%lT, state.getPotentialEnergy())
    positions = state.getPositions()

  # Save checkpoint after temperature ramp
  simulation.context.setPositions(positions)
  state = simulation.context.getState(getPositions=True, getVelocities=True, \
    enforcePeriodicBox=True)
  checkpointReporter.report(simulation, state)
  restartReporter.report(simulation, state)
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
simulation.reporters.append(restartReporter)
simulation.step(steps)

# Store checkpoint and restart
state = simulation.context.getState(getPositions=True, getVelocities=True, \
  enforcePeriodicBox=True)
checkpointReporter.report(simulation, state)
restartReporter.report(simulation, state)
