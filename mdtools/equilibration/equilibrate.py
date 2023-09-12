from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj
from mdtraj.reporters import HDF5Reporter

def equilibrate(mdsystem, simtime=1.*nanoseconds, temperature=300*kelvin, posre=True, reportLargeForceThreshold=-1):
    """
    Minimizes and equilibrate an MDSystem object. If position restraints
    are applied, it will taper the restraints over the course of the 
    simulation. This method assumes that MDSystem.buildSimulation() has
    already been called.

    Half of the simulation time is used for tapering restraints, and the
    other half is simulated without position restraints.

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
    reportLargeForceThreshold : int
            If <= 0, will not report; otherwise, print a list of
            all atoms with net forces on them exceeding the
            threshold in magnitude
    """

    # Minimize system
    mdsystem.minimize()

    # Initialize velocities
    mdsystem.simulation.context.setVelocitiesToTemperature(temperature)

    # Taper restraints if they're applied
    if posre:
        splits = 41
        for i in range(splits):
            k = max((5.0 - (0.25*i)), 0)
            mdsystem.simulation.context.setParameter('k', k*kilocalories_per_mole/angstroms**2)
            mdsystem.simulate(simtime/splits, reportLargeForceThreshold=reportLargeForceThreshold)
    else:
        mdsystem.simulate(simtime, reportLargeForceThreshold=reportLargeForceThreshold)
            
    return mdsystem
