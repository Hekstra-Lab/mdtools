#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from mdtools.prep import supercell

def buildUnitCell(filename, spacegroup=None, cutoff=3.2*angstroms):
    """
    Build up unitcell for simulation
    
    Parameters
    ----------
    spacegroup : str
        Crystallographic space group for building up the unit cell
    cutoff : distance
        Distance threshold for considering waters to be clashing
    """
    pdb = PDBFile(filename)
    m = Modeller(pdb.topology, pdb.positions)
    unitcell = supercell.buildUnitCell(m, spacegroup)
    unitcell = supercell._imageUnitCell(unitcell)
    ff = ForceField('amber99sbildn.xml', 'tip3p.xml')
    preppedUC = supercell._solvateUnitCell(unitcell, ff, 
                                           distThreshold=cutoff, 
                                           ionicStrength=0.1*molar)
    return preppedUC

if __name__ == "__main__":

    preppedUC = buildUnitCell("pdb/3ZEK_prepped.pdb", "P 43 21 2", cutoff=2.8*angstroms)

    # Write out PDB of system
    with open("3ZEK_unitcell.pdb", "w") as pdbout:
        PDBFile.writeFile(preppedUC.getTopology(), preppedUC.getPositions(), pdbout)

