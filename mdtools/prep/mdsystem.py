"""
mdsystem.py: Extends Modeller class from OpenMM with useful methods

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
import numpy as np
from mdtools.equilibration import equilibrate

class MDSystem(Modeller):
    """
    MDSystem extends the Modeller class, providing useful methods for 
    preparing molecular systems for MD simulations in OpenMM.
    """
    def __init__(self, topology, positions, forcefield):
        Modeller.__init__(self, topology, positions)
        self.forcefield = forcefield

    def _toMDTrajTopology(self):
        """
        Returns a MDTraj Topology object from the OpenMM Topology of 
        this system. Importantly, it also ensures that the resSeq
        attribute is set. This is necessary to ensure that PDB-based
        residue numbering is not lost.

        Returns
        -------
        mdtrajtop : MDTraj.Topology
            MDTraj Topology object representing the system
        """
        mdtrajtop = mdtraj.Topology.from_openmm(self.getTopology())
        for r1, r2 in zip(self.topology.residues(), mdtrajtop.residues):
            r2.resSeq = int(r1.id)
        return mdtrajtop

    def findMolecules(self):
        """
        Identify molecules based on bonded subsets of atoms.

        Returns
        -------
        molecules : list of sets of atoms
            Each entry represents one molecule, and is the set of all 
            Atoms in that molecule
        """
        mdtrajtop = self._toMDTrajTopology()
        return mdtrajtop.find_molecules()

    def select(self, selection):
        """
        Select the atoms in the system that match the provided selection
        string. Uses the MDTraj syntax to define atom selections.

        Returns
        -------
        indices : np.ndarray
            Array of the indices of atoms matching the selection
        """
        mdtrajtop = self._toMDTrajTopology()
        return mdtrajtop.select(selection)
        
    def buildSimulation(self, temperature=300*kelvin, ensemble="NPT", posre=False,
                        nonbondedMethod=PME, nonbondedCutoff=1.*nanometer):
        """
        Build a simulation context from the system. The simulation is
        then available as an attribute.
        """
        # Build system
        system = self.forcefield.createSystem(self.topology, nonbondedMethod=nonbondedMethod, 
                                              nonbondedCutoff=nonbondedCutoff, 
                                              constraints=HBonds)

        # Setup MD simulation
        dt = 0.002*picoseconds
        integrator = LangevinIntegrator(temperature, 1/picosecond, dt)

        # Add position restraints that can be tapered off during simulation
        if posre:
            force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            force.addGlobalParameter("k", 5.0*kilocalories_per_mole/angstroms**2)
            force.addPerParticleParameter("x0")
            force.addPerParticleParameter("y0")
            force.addPerParticleParameter("z0")
            for i in self.select("not water and not (element Na or element Cl) and not element H"):
                force.addParticle(int(i), positions[i].value_in_unit(nanometers))
            system.addForce(force)
        
        # Setup barostat for NPT ensemble
        if ensemble == "NPT":
            barostat   = MonteCarloBarostat(1.0*bar, temperature, 25)
            system.addForce(barostat)

        # Add simulation
        self.simulation = Simulation(self.topology, system, integrator)

        # Initialize particle positions
        self.simulation.context.setPositions(self.positions)

        return self

    def minimize(self):
        """
        Minimize the system using the simulation context. If a simulation
        context has not been built, an attribute error is raised.
        """
        self.simulation.minimizeEnergy()

        # Update positions
        state = self.simulation.context.getState(getPositions=True)
        self.positions = state.getPositions()

        return self

    def _time2steps(self, time):
        """
        Compute the number of steps corresponding to a given chemical time
        """
        chemtime = time.in_units_of(picoseconds)
        dt = self.simulation.integrator.getStepSize()
        return int(np.ceil(chemtime / dt))
        
    def simulate(self, n):
        """
        Simulate the system for the given number of steps. If n is a 
        simtk.Unit of time, the number of steps are chosen to simulate
        for the indicated chemical time

        Parameters
        ----------
        n : int or simtk.unit
            Number of steps or chemical time of simulation
        """
        if isinstance(n, int):
            self.simulation.step(n)
        else:
            self.simulation.step(self._time2steps(n))

        # Update positions
        state = self.simulation.context.getState(getPositions=True)
        self.positions = state.getPositions()

        return self

    def equilibrate(self, simtime=1.*nanoseconds, temperature=300*kelvin, posre=True):
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
        return equilibrate.equilibrate(self, simtime, temperature, posre)
