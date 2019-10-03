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
from scipy import stats
import pandas as pd
import itertools

def squeeze(mdsystem, tolerance=0.003, maxIterations=10):
    """
    Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors
    """
    # Target periodic box vectors
    targetv = mdsystem.topology.getPeriodicBoxVectors()
    target = np.array(targetv.value_in_unit(nanometers))
    targetvol = np.linalg.det(target)
    
    # Loop until converged
    iteration = 1
    while iteration <= maxIterations:

        # Calmdown to relieve bad steric clashes
        mdsystem.calmdown(posre=True)
        
        # Build simulation
        mdsystem.buildSimulation(ensemble="NPT", posre=True, filePrefix=f"iter{iteration:02d}",
                                 saveTrajectory=True, saveStateData=True,
                                 trajInterval=1000, stateDataInterval=1000)

        # Save initial positions
        state = mdsystem.simulation.context.getState(getPositions=True)
        startingPositions = np.array(state.getPositions().value_in_unit(nanometers))

        # Equilibration run, tapering off position restraints
        mdsystem.equilibrate(simtime=2.0*nanoseconds, posre=True)

        # Assess convergence of unit cell volume
        simtime = 250
        converged = False
        while True:

            # Compute mean and standard error of unit cell volume
            df = pd.read_csv(f"iter{iteration:02d}.csv")
            vol   = df["Box Volume (nm^3)"].iloc[-simtime:].mean() 
            stderr = stats.sem(df["Box Volume (nm^3)"].iloc[-simtime:])
            
            # Check convergence criteria
            percent_diff1 = (targetvol - (vol+stderr)) / targetvol
            percent_diff2 = (targetvol - (vol-stderr)) / targetvol
            change = targetvol - vol

            # Case 1: Both error bounds are on one side of targetvol
            if (((percent_diff1 < 0) and (percent_diff2 < 0)) or
                ((percent_diff1 > 0) and (percnet_diff2 > 0))):
                converged = False
                break
            
            # Case 2: stderr is too high 
            elif ((stderr/vol) > tolerance):
                mdsystem.simulate(1.0*nanoseconds)
                simtime += 500
                continue

            # Case 3: Too much simulation time -- probably not correct
            elif simtime > 5000:
                converged = False
                break
            
            # Case 4: simulation cell has converged within error margins
            elif (np.abs(percent_diff1) < tolerance) and (np.abs(percent_diff2) < tolerance):
                converged = True
                break
            
            else:
                converged = False
                break

        print(f"Percent Change: {change/targetvol} +/- {stderr/targetvol} ")

        # Close open files
        for reporter in mdsystem.simulation.reporters:
            try:
                reporter.close()
            except:
                continue

        # Check squeeze run convergence criteria
        if converged:
            break

        # Revert positions and box vectors
        mdsystem.positions = startingPositions*nanometers
        mdsystem.topology.setPeriodicBoxVectors(targetv)

        # Add or delete waters
        numWaters = np.floor(np.abs(change)/0.1)
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
    randwaters = np.random.choice(chain.n_residues, min(numWaters, chain.n_residues), replace=False)
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

    # Select n random waters to delete
    randwaters = np.random.choice(chain.n_residues, min(numWaters, chain.n_residues), replace=False)
    residues   = [ chain.residue(i) for i in randwaters ]
    atoms = itertools.chain(*[ r.atoms for r in residues ])
    atomindices  = [ a.index for a in atoms ]
    atoms = [ a for a in mdsystem.topology.atoms() if a.index in atomindices ]
    mdsystem.delete(atoms)

    return
