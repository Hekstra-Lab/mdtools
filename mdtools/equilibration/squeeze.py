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
import numpy as np

def squeeze(mdsystem, tolerance=0.003):
    """
    Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors
    """
    # Target periodic box vectors
    target = mdsystem.topology.getPeriodicBoxVectors()
    target = np.array(target.value_in_unit(nanometers))

    # Short pre-equilibration run
    mdsystem.buildSimulation(ensemble="NPT", posre=True, filePrefix="pre_equilibration",
                             saveTrajectory=True, saveStateData=True)
    mdsystem.equilibrate(simtime=100*picoseconds, posre=False)
    
    # Loop until converged
    iteration = 1
    while True:
    
        # Drop checkpoint
        mdsystem.buildSimulation(ensemble="NPT", posre=True, filePrefix=f"iter{iteration}",
                             saveTrajectory=True, saveStateData=True)
        checkpt = mdsystem.simulation.context.createCheckpoint()
    
        # Equilibration run, tapering off position restraints
        mdsystem.equilibrate(simtime=1.*nanoseconds, posre=True)
        
        # Determine change in periodic box vectors
        state   = mdsystem.simulation.context.getState() 
        current = state.getPeriodicBoxVectors()
        current = np.array(current.value_in_unit(nanometers))
        change = (target - current).max()

        # Convergence criteria
        if change < tolerance:
            break
        
        # Load checkpoint
        mdsystem.simulation.context.loadCheckpoint(checkpt)

        # Add or delete waters
        numWaters = np.floor(np.abs(change)/0.00147)
        if change > 0.0:
            for i in range(numWaters):
                duplicateWater(mdsystem)
        else:
            for i in range(numWaters):
                deleteWater(mdsystem)

        # Increment iteration
        iteration += 1
                
    return

def duplicateWater(mdsystem):
    """
    Duplicate a random water
    """
    # Select waters
    waters = mdsystem.select("water")
    mdtrajtop    = mdtraj.Topology.from_openmm(mdsystem.topology)
    mdtrajwaters = mdtrajtop.subset(waters)

    # Select a random water to duplicate
    residue      = mdtrajwaters.residue(np.random.randint(len(waters)/3))
    atomindices  = [ a.index for a in residue.atoms ]
    newtop    = mdtrajwaters.subset(atomindices)
    newtop = newtop.to_openmm()
    pos = np.array(mdsystem.positions.value_in_unit(nanometers))

    # Jitter the position to avoid singularity and add to system
    jitter = (0.5-np.random.rand(3))*0.1
    mdsystem.add(newtop, (pos[atomindices] + jitter)*nanometers)

    return

def deleteWater(mdsystem):
    """
    Delete a random water
    """
    # Select waters
    waters = mdsystem.select("water")
    mdtrajtop    = mdtraj.Topology.from_openmm(mdsystem.topology)
    mdtrajwaters = mdtrajtop.subset(waters)

    # Select a random water to delete
    residue      = mdtrajwaters.residue(np.random.randint(len(waters)/3))
    atomindices  = [ a.index for a in residue.atoms ]
    atoms = [ a for a in mdsystem.topology.atoms() if a.index in atomindices ]
    mdsystem.delete(atoms)

    return
