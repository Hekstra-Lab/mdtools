"""
latticemdsystem.py: Provides LatticeMDSystem class for preparing MD
                    simulations of crystal lattices

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
from .mdsystem import MDSystem

class LatticeMDSystem(MDSystem):
    """
    LatticeMDSystem provides methods for preparing MD simulations of 
    crystal lattices in OpenMM.
    """
    def __init__(self, topology, positions, forcefield, spacegroup):
        MDSystem.__init__(self, topology, positions, forcefield)
        self.spacegroup = spacegroup
