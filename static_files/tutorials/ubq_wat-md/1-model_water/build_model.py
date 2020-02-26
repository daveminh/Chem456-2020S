from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

pdb = app.PDBFile('../0-propka/1ubq.pqr')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
modeller.addSolvent(forcefield, model='tip3p', padding=1*unit.nanometers, \
  positiveIon='Na+', negativeIon='Cl-', \
  ionicStrength=unit.Quantity(value=0.150, unit=unit.molar), \
  neutralize=True)
pdb = modeller
app.PDBFile.writeFile(pdb.topology, pdb.positions, open('1ubq.pdb','w'))
