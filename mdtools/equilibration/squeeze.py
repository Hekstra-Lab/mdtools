"""
squeeze.py: Squeeze run to titrate the number of waters in lattice

Author: Jack Greisman <greisman@g.harvard.edu>
        Ziyuan Zhao   <ziyuanzhao@fas.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.1"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
import numpy as np
from scipy import stats
import pandas as pd
import itertools

def _iter_duplicate_waters(mdsystem, n_waters, dn = 1000):
    if n_waters == 0:
        return
    if dn <= 0:
        duplicateWaters(mdsystem, n_waters)
        return
    while True:
        if dn <= n_waters:
            n_waters -= dn
            duplicateWaters(mdsystem, dn)
            mdsystem.calmdown(posre=True)
        else:
            duplicateWaters(mdsystem, n_waters)
            break

def _iter_delete_waters(mdsystem, n_waters, dn = 1000):
    if n_waters == 0:
        return
    if dn <= 0:
        deleteWaters(mdsystem, n_waters)
        return
    while True:
        if dn <= n_waters:
            n_waters -= dn
            deleteWaters(mdsystem, dn)
            mdsystem.calmdown(posre=True)
        else:
            deleteWaters(mdsystem, n_waters)
            break

"""Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors.


    Args:
        mdsystem (mdtools.MDSystem): System to perform the squeeze run on
        tolerance (float, optional): Threshold for volume change after squeeze iteration for determining whether the system has reached convergence. Defaults to 0.003.
        maxIterations (int, optional): Maximum iterations of titration to perform. Defaults to 10.
        maxSimtime (simtk.Unit, optional): Maximum total simulation time allowed for the squeeze run. Defaults to 10*nanoseconds.
        initial_water_perturb (int, optional): Initial number of water molecules to add to the system. Defaults to 0.
        dn (int, optional): Number of water molecules to add each time. If set to 0, will add all water molecules in one step. Defaults to 0.
        dt (simtk.Unit, optional): Timestep of the simulation used for squeezing the system. Defaults to 0.002*picoseconds.
        targetvol (int, optional): Target volume. If set to 0, will use the initial volume of the input system. Defaults to 0.
    """

def squeeze(mdsystem, tolerance=0.003, maxIterations=10, maxSimtime=10*nanoseconds, initial_water_perturb = 0, dn = 0, dt=0.002*picoseconds, targetvol=0):
    """Squeeze run to titrate the number of waters in a lattice MD system
    in order to maintain desired periodic box vectors.

    Parameters
    ----------
    mdsystem : mdtools.MDSystem
        System to perform the squeeze run on
    tolerance : float, optional
        Threshold for volume change after squeeze iteration for determining
        whether the system has reached convergence, by default 0.003
    maxIterations : int, optional
        Maximum iterations of titration to perform, by default 10
    maxSimtime : simtk.Unit, optional
        Maximum total simulation time allowed for the squeeze run, 
        by default 10*nanoseconds
    initial_water_perturb : int, optional
        Initial number of water molecules to add to the system, 
        If set to 0, will add all water molecules in one step, by default 0
    dn : int, optional
        Number of water molecules to add each time, by default 0
    dt : simtk.Unit, optional
        Timestep of the simulation used for squeezing the system,
        by default 0.002*picoseconds
    targetvol : int, optional
        Target volume. If set to 0, will use the initial volume of the 
        input system, by default 0
    """

    # Target periodic box vectors
    targetv = mdsystem.topology.getPeriodicBoxVectors()

    if targetvol==0:
        targetv = mdsystem.topology.getPeriodicBoxVectors()
        target = np.array(targetv.value_in_unit(nanometers))
        targetvol = np.linalg.det(target)

    # Initial addition or removal of water if required
    if initial_water_perturb > 0:
        _iter_duplicate_waters(mdsystem, initial_water_perturb, dn)
        print(f"Initial perturbation +{initial_water_perturb} waters", flush=True)
    elif initial_water_perturb < 0:
        _iter_delete_waters(mdsystem, -initial_water_perturb, dn)
        print(f"Initial perturbation -{-initial_water_perturb} waters", flush=True)

    # Loop until converged
    iteration = 1
    while iteration <= maxIterations:

        # Calmdown to relieve bad steric clashes
        mdsystem.calmdown(posre=True)
        
        # Build simulation
        mdsystem.buildSimulation(ensemble="NPT", posre=True, filePrefix=f"iter{iteration:02d}",
                                 saveTrajectory=True, saveStateData=True,
                                 trajInterval=1000, stateDataInterval=1000,
				                 dt=dt)

        # Save initial positions
        state = mdsystem.simulation.context.getState(getPositions=True)
        startingPositions = np.array(state.getPositions().value_in_unit(nanometers))

        # Equilibration run, tapering off position restraints
        mdsystem.equilibrate(simtime=2.0*nanoseconds, posre=True)

        # Case 1: Too far from convergence, adjust number of waters
        vol, sem, std = _computeVolumeStats(f"iter{iteration:02d}.csv", 250)
        change = targetvol - vol
        if (change / targetvol) > 3*tolerance:
            converged = False
            print("Case 1: Too far from convergence", flush=True)

        # Otherwise, simulate more, and assess convergence
        else:
            mdsystem.simulate(maxSimtime - 2.0*nanoseconds)
            steps = int(np.ceil(maxSimtime/dt)) - 750
            vol, sem, std = _computeVolumeStats(f"iter{iteration:02d}.csv", steps)
            change = targetvol - vol
            
            # Check convergence criteria
            percent_diff1 = (targetvol - (vol+sem)) / targetvol
            percent_diff2 = (targetvol - (vol-sem)) / targetvol

            # Case 2: simulation cell has converged within error margins
            if ((np.abs(percent_diff1) < tolerance) and (np.abs(percent_diff2) < tolerance)):
                converged = True
                print("Case 2: converged", flush=True)

            else:
                converged = False
                print("Case 3: simulated more; not converged", flush=True)

        print(f"Percent Change: {change/targetvol} +/- {sem/targetvol}", flush=True)

        # Close open files
        for reporter in mdsystem.simulation.reporters:
            try:
                reporter.close()
            except:
                continue

        # Check squeeze run convergence criteria
        if converged:
            mdsystem.saveCheckpoint("squeezed.chkpt")
            break

        # Revert positions and box vectors
        mdsystem.positions = startingPositions*nanometers
        mdsystem.topology.setPeriodicBoxVectors(targetv)
        print("Reverted atom positions and box vectors")

        # Add or delete waters
        numWaters = np.floor(np.abs(change)/0.05)
        if change > 0.0:
            iter_duplicate_waters(mdsystem, int(numWaters), dn)
            print(f"+{numWaters} waters", flush=True)
        else:
            iter_delete_waters(mdsystem, int(numWaters), dn)
            print(f"-{numWaters} waters", flush=True)

        # Increment iteration
        iteration += 1

    return

def _computeVolumeStats(csvfile, simtime):
    """
    Determine mean, sem, and std of box volume. 

    Parameters
    ----------
    csvfile : str
        CSV file from OpenMM StateDataReporter
    simtime : int
        Number of frames from simulation end to consider for 
        box volume statistics

    Returns
    -------
    (vol, sem, std) : tuple
         Mean, standard error, and standard deviation of the box volume
    """
    df = pd.read_csv(csvfile)
    vol = df["Box Volume (nm^3)"].iloc[-simtime:].mean()
    sem = stats.sem(df["Box Volume (nm^3)"].iloc[-simtime:])
    std = df["Box Volume (nm^3)"].iloc[-simtime:].std()
    return vol, sem, std
    
def duplicateWaters(mdsystem, numWaters):
    """Duplicate a random water.

    Parameters
    ----------
    mdsystem : mdtools.MDSystem
        The system to duplicate water molecules from
    numWaters : int
        Number of individual water molecules to duplicate
    """

    # Find chain with most waters
    mdtrajtop    = mdsystem._toMDTrajTopology()
    watercounts = []
    for c in mdtrajtop.chains:
        n = sum([ 1 for r in c.residues if r.name == "HOH" ])
        watercounts.append(n)
    chain = mdtrajtop.chain(np.argmax(watercounts))

    # Select n random waters to duplicate
    water_idxs = [i for i in range(chain.n_residues) if chain.residue(i).name == "HOH"] # this avoids grabbing nonwater residues
    randwaters = np.random.choice(water_idxs, min(numWaters, len(water_idxs)), replace=False)
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
    """Delete a random water.

    Parameters
    ----------
    mdsystem : mdtools.MDSystem
        The system to delete water molecules from
    numWaters : int
        Number of individual water molecules to delete
    """

    # Find chain with most waters
    mdtrajtop    = mdsystem._toMDTrajTopology()
    watercounts = []
    for c in mdtrajtop.chains:
        n = sum([ 1 for r in c.residues if r.name == "HOH" ])
        watercounts.append(n)
    chain = mdtrajtop.chain(np.argmax(watercounts))

    # Select n random waters to delete
    water_idxs = [i for i in range(chain.n_residues) if chain.residue(i).name == "HOH"] # this avoids grabbing nonwater residues
    randwaters = np.random.choice(water_idxs, min(numWaters, len(water_idxs)), replace=False)
    residues   = [ chain.residue(i) for i in randwaters ]
    atoms = itertools.chain(*[ r.atoms for r in residues ])
    atomindices  = [ a.index for a in atoms ]
    atoms = [ a for a in mdsystem.topology.atoms() if a.index in atomindices ]
    mdsystem.delete(atoms)

    return
