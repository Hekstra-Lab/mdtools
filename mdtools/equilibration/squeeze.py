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

    mdsystem.buildSimulation(ensemble="NPT", posre=True)
    
    # Short pre-equilibration run

    # Drop checkpoint

    #----------------------------------------------------------------------#
    # LOOP until converged

    # Load checkpoint
    
    # Equilibration run, tapering off position restraints

    # Determine change in periodic box vectors

    # Duplicate waters or delete waters
    #----------------------------------------------------------------------#

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
    newtop    = mdtrajwaterssubset(atomindices)
    newtop = newtop.to_openmm()
    pos = np.array(mdsystem.positions.value_in_unit(nanometers))

    # Jitter the position to avoid singularity and add to system
    jitter = (0.5-np.random.rand(3))*0.1
    mdsystem.add(newtop, (pos[atomindices] + jitter)*nanometers)

    return
