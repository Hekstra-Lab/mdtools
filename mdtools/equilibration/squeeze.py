"""
squeeze.py: Squeeze run to titrate the number of waters in lattice

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj

def squeeze(mdsystem):
    """
    Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors
    """

    return
