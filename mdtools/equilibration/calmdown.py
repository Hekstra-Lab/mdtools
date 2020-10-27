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
from itertools import product
from scipy.spatial.distance import cdist

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
    # Build simulation with Brownian dynamics
    for i in range(10):
        mdsystem.buildSimulation(ensemble="NVT", posre=posre)
        problemPairs = _identifyProblemPairs(mdsystem)
        mdsystem.buildSimulation(integrator=BrownianIntegrator, dt=0.001*attoseconds,
                                 ensemble="NVT", exceptions=problemPairs, posre=posre,
                                 constraints=None, rigidWater=False)
        mdsystem.simulate(10000)
    mdsystem.minimize()

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
    state = mdsystem.simulation.context.getState(getPositions=True,
                                                 getForces=True)
    netforces = np.linalg.norm(state.getForces(asNumpy=True), axis=1)
    indices = np.where(np.isnan(netforces) | (netforces > 5e3))[0]

    # Return list of force overflow atoms that are also < 5A from eachother
    # using periodic distance
    positions = state.getPositions(asNumpy=True)[indices]
    positions = np.array(positions.value_in_unit(nanometers))
    pos = positions[np.newaxis, :, :]
    vectors = state.getPeriodicBoxVectors()
    dims = internal.unitcell.computeLengthsAndAngles(vectors)
    neighborhood = np.repeat(pos, 27, axis=0)
    cells = [-1, 0, 1]
    transmat = np.array([ uc for uc in product(cells, cells, cells) ])*dims[:3]
    neighborhood += transmat.reshape(27, 1, 3)
    dists = np.array([ cdist(positions, neighbor) for neighbor in neighborhood ])
    minimage_dist = dists.min(axis=0)
    pairs = [ (indices[i], indices[j]) for  i, j in zip(*np.where(minimage_dist < 0.5)) if i != j ]
    
    return pairs
    
