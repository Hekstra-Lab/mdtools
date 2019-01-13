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
        
    def addSolvent(self, **kwargs):
        """
        Add solvent to the system to fill a rectangular box.

        Parameters
        ----------
        model : str='tip3p'
            the water model to use.  Supported values are 'tip3p', 'spce', 'tip4pew', and 'tip5p'.
        boxSize : Vec3=None
            the size of the box to fill with water
        boxVectors : tuple of Vec3=None
            the vectors defining the periodic box to fill with water
        padding : distance=None
            the padding distance to use
        numAdded : int=None
            the total number of molecules (waters and ions) to add
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        """
        Modeller.addSolvent(self, self.forcefield, **kwargs)
        
    def addMembrane(self, **kwargs):
        """
        Add a lipid membrane to the system.

        Parameters
        ----------
        lipidType : string or object
            the type of lipid to use.  Supported string values are 'POPC' and 'POPE'.  For other types
            of lipids, provide a PDBFile or PDBxFile object (or any other object with "topology" and
            "positions" fields) containing a membrane patch.
        membraneCenterZ: distance=0*nanometer
            the position along the Z axis of the center of the membrane
        minimumPadding : distance=1*nanometer
            the padding distance to use
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        """
        Modeller.addMembrane(self, self.forcefield, **kwargs)
