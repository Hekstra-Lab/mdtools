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
import itertools

def calmdown(mdsystem, posre=True):
    """
    Aggressive relaxation of a molecular system for MD system. 

    Protocol:
        1) Clashes/overlapping positions are identified, and their 
           nonbonded interactions are excluded to prevent force 
           overflows. 
        2) Brownian dynamics is then used to gently equilibrate the 
           system using a very short time step and without constraints
           on hydrogens or water
        3) The exceptions are then removed and Brownian dynamics is 
           repeated with a slightly longer timestep, hydrogen constraints,
           and rigid waters
        4) Finally, Langevin dynamics is simulated with a 2 fs timestep
           to ensure the simulation can be simulated.

    Parameters
    ----------
    mdsystem : MDSystem
        Molecular system to be relaxed
    posre : bool
        If true, position restraints are applied to non-water, non-ion
        heavy atoms in the system.
    """

    # Identify problematic atom pairs
    mdsystem.buildSimulation(ensemble="NVT", posre=posre)
    problemPairs = _identifyProblemPairs(mdsystem)

    # Build simulation with Brownian dynamics and exceptions
    mdsystem.buildSimulation(integrator=BrownianIntegrator, dt=0.001*attoseconds,
                             ensemble="NVT", exceptions=problemPairs, posre=posre,
                             constraints=None, rigidWater=False)
    mdsystem.simulate(10000)
    
    # Build simulation with Brownian dynamics
    mdsystem.buildSimulation(integrator=BrownianIntegrator, dt=0.001*femtoseconds,
                             ensemble="NVT", posre=posre)
    mdsystem.minimize()
    mdsystem.simulate(10*femtoseconds)

    # Build simulation with Langevin dynamics
    mdsystem.buildSimulation(integrator=LangevinIntegrator, dt=2*femtoseconds,
                             ensemble="NVT", posre=posre)
    mdsystem.minimize()
    mdsystem.simulate(1*picoseconds)
    
    return

def _identifyProblemPairs(mdsystem):
    """
    Identify problematic pairs of atoms based on overflows in force
    calculations

    Parameters
    ----------
    mdsystem : MDSystem
        Molecular system to inspect for problematic pairs of atoms
    
    Returns
    -------
    problemPairs : list of tuples
        List of (atom1, atom2) pairs that should be excluded from
        nonbonded force calculations
    """
    # Identify atoms with force overflows or large forces
    state = mdsystem.simulation.context.getState(getForces=True)
    netforces = np.linalg.norm(state.getForces(asNumpy=True), axis=1)
    indices = np.where(np.isnan(netforces))[0]
    
    return itertools.combinations(indices, 2)
    
