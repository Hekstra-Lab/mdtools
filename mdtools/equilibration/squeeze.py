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
import itertools

def squeeze(mdsystem, tolerance=0.003, maxIterations=10):
    """
    Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors
    """
    # Target periodic box vectors
    targetv = mdsystem.topology.getPeriodicBoxVectors()
    target = np.array(targetv.value_in_unit(nanometers))
    
    # Loop until converged
    iteration = 1
    while iteration <= maxIterations:

        # Build simulation
        mdsystem.buildSimulation(ensemble="NPT", posre=True, filePrefix=f"iter{iteration}",
                                 saveTrajectory=True, saveStateData=True,
                                 trajInterval=1000, stateDataInterval=1000)

        # Save initial positions
        state = mdsystem.simulation.context.getState(getPositions=True)
        startingPositions = np.array(state.getPositions().value_in_unit(nanometers))

        # Equilibration run, tapering off position restraints
        mdsystem.equilibrate(simtime=2.0*nanoseconds, posre=True)

        # Close open files
        for reporter in mdsystem.simulation.reporters:
            try:
                reporter.close()
            except:
                continue

        # Determine change in periodic box vectors
        traj = mdtraj.load(f"iter{iteration}.h5")
        current = np.mean(traj.unitcell_lengths[-100:], axis=0)
        percent_change = (np.abs(np.diag(target) - current)/np.diag(target)).max()
        change = (np.diag(target) - current).max()
        print(f"Percent Change: {percent_change}")

        # Convergence criteria
        if percent_change < tolerance:
            break

        # Revert positions and box vectors
        mdsystem.positions = startingPositions*nanometers
        mdsystem.topology.setPeriodicBoxVectors(targetv)

        # Add or delete waters
        numWaters = np.floor(np.abs(change)/0.001)
        if change > 0.0:
            duplicateWaters(mdsystem, int(numWaters))
        else:
            deleteWaters(mdsystem, int(numWaters))

        # Increment iteration
        iteration += 1

    return

def duplicateWaters(mdsystem, numWaters):
    """
    Duplicate a random water
    """
    # Find chain with most waters
    mdtrajtop    = mdsystem._toMDTrajTopology()
    watercounts = []
    for c in mdtrajtop.chains:
        n = sum([ 1 for r in c.residues if r.name == "HOH" ])
        watercounts.append(n)
    chain = mdtrajtop.chain(np.argmax(watercounts))

    # Select n random waters to duplicate
    randwaters = np.random.choice(chain.n_residues, numWaters, replace=False)
    residues   = [ chain.residue(i) for i in randwaters ]
    atoms = itertools.chain(*[ r.atoms for r in residues ])
    atomindices  = [ a.index for a in atoms ]
    newtop    = mdtrajtop.subset(atomindices)
    newtop = newtop.to_openmm()
    pos = np.array(mdsystem.positions.value_in_unit(nanometers))

    # Jitter the position to avoid singularity and add to system
    jitter = (0.5-np.random.rand(len(atomindices), 3))*0.2
    mdsystem.add(newtop, (pos[atomindices] + jitter)*nanometers)

    return

def deleteWaters(mdsystem, numWaters):
    """
    Delete a random water
    """
    # Find chain with most waters
    mdtrajtop    = mdsystem._toMDTrajTopology()
    watercounts = []
    for c in mdtrajtop.chains:
        n = sum([ 1 for r in c.residues if r.name == "HOH" ])
        watercounts.append(n)
    chain = mdtrajtop.chain(np.argmax(watercounts))

    # Select a random water to delete
    randwaters = np.random.choice(chain.n_residues, numWaters, replace=False)
    residues   = [ chain.residue(i) for i in randwaters ]
    atoms = itertools.chain(*[ r.atoms for r in residues ])
    atomindices  = [ a.index for a in atoms ]
    atoms = [ a for a in mdsystem.topology.atoms() if a.index in atomindices ]
    mdsystem.delete(atoms)

    return
