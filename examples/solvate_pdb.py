"""
solvate_pdb.py
Jack Greisman

Example script to demonstrate how to use OpenMM to solvate a prepped PDB
file for simulation
"""

from simtk.openmm.app import PDBFile
from simtk.openmm.app import Modeller
from simtk.openmm.app import ForceField
from simtk.unit import *

# Load PDB and instantiate a Modeller object
pdb = PDBFile("pdb/3ZEK_prepped.pdb")
m = Modeller(pdb.getTopology(), pdb.getPositions())

# Solvate a system
ff = ForceField("amber99sbildn.xml", "tip3p.xml")
m.addSolvent(ff, padding=15*angstroms, ionicStrength=0.3*molar, neutralize=True)

# Write out PDB of system
with open("3ZEK_solvated.pdb", "w") as pdbout:
    PDBFile.writeFile(m.getTopology(), m.getPositions(), pdbout)
