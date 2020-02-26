from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import os

# Input files

system_name = '1ubq'
pdb = app.PDBFile('../1-model_water/'+system_name+'.pdb')
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
while os.path.isfile('%s%03d.dcd'%(system_name,run_number)):
  run_number += 1
run_name = '%s%03d'%(system_name,run_number)
print('Data will be stored with the prefix '+run_name)

steps = int(100*unit.nanoseconds/dt)
print('Simulation will run for %d steps\n'%steps)

if os.path.isdir('/Users/dminh/'):
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
checkpointReporter = app.CheckpointReporter(system_name+'.chk', 5000)

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

if not os.path.isfile(system_name+'.chk'):
  print('Performing energy minimization...')
  simulation.minimizeEnergy()
  simulation.currentStep = 0
else:
  print('Loading from checkpoint')
  with open(system_name+'.chk', 'rb') as f:
    simulation.context.loadCheckpoint(f.read())

# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.step(steps-simulation.currentStep)
