from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj
from mdtraj.reporters import HDF5Reporter

def equilibrate(mdsystem, simtime=1.*nanoseconds, temperature=300*kelvin, posre=True):
    """
    Minimizes and equilibrate an MDSystem object. If position restraints
    are applied, it will taper the restraints over the course of the 
    simulation. This method assumes that MDSystem.buildSimulation() has
    already been called.

    Parameters
    ----------
    mdsystem : MDSystem
        MDSystem object to equilibrate
    simtime : simtk.unit
        Total simulation time to use for equilibration
    temperature : OpenMM.unit.Quantity(unit=kelvin)
        Temperature to use to initialize velocities
    posre : bool
        If True, position restraints have been applied to simulation object
    """

    # Minimize system
    mdsystem.minimize()

    # Initialize velocities
    mdsystem.simulation.context.setVelocitiesToTemperature(temperature)

    # Taper restraints if they're applied
    if posre:
        for i in range(11):
            k = max((5.0 - (0.5*i)), 0)
            mdsystem.simulation.context.setParameter('k', k*kilocalories_per_mole/angstroms**2)
            mdsystem.simulate(simtime/11)
    else:
        mdsystem.simulate(simtime)
            
    return mdsystem
