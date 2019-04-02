"""
calmdown.py: Use Brownian dynamics and minimization to relax a 
             molecular system for MD simulations.

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
import numpy as np

def calmdown(mdsystem, posre=True):
    """
    Aggressive relaxation of a molecular system for MD system. 

    Protocol:
        1) Clashes/overlapping positions are identified, and their 
           nonbonded interactions are excluded to prevent force 
           overflows. 
        2) Brownian dynamics is then used to gently equilibrate the 
           system. 
        3) The exclusions are then removed and Brownian dynamics is 
           repeated.
        4) Local energy minimizer is used to relax system.

    Parameters
    ----------
    mdsystem : MDSystem
        Molecular system to be relaxed
    posre : bool
        If true, position restraints are applied to non-water, non-ion
        heavy atoms in the system.
    """

    
    return
