"""
solvatedmdsystem: Provides SolvatedMDSystem class for preparing MD
                  simulations of solvated systems

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
from .mdsystem import MDSystem

class SolvatedMDSystem(MDSystem):
    """
    SolvatedMDSystem provides methods for preparing MD simulations of 
    solvated systems
    """
    def __init__(self, topology, positions, forcefield):
        MDSystem.__init__(self, topology, positions, forcefield)
