"""
mdsystem.py: Extends Modeller class from OpenMM with useful methods

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj

class MDSystem(Modeller):
    """
    MDSystem extends the Modeller class, providing useful methods for 
    preparing molecular systems for MD simulations in OpenMM.
    """
    def __init__(self, topology, positions, forcefield):
        Modeller.__init__(self, topology, positions)
        self.forcefield = forcefield

    def findMolecules(self):
        """
        Identify molecules based on bonded subsets of atoms.

        Returns
        -------
        molecules : list of sets of atoms
            Each entry represents one molecule, and is the set of all 
            Atoms in that molecule
        """
        mdtrajtop = mdtraj.Topology.from_openmm(self.getTopology())
        return mdtrajtop.find_molecules()

    def select(self, selection):
        """
        Select the atoms in the system that match the provided selection
        string. Uses the MDTraj syntax to define atom selections.

        Returns
        -------
        indices : np.ndarray
            Array of the indices of atoms matching the selection
        """
        mdtrajtop = mdtraj.Topology.from_openmm(self.getTopology())
        return mdtrajtop.select(selection)
        
