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
from mdtraj.reporters import HDF5Reporter
import numpy as np
from mdtools.equilibration import equilibrate, calmdown

class MDSystem(Modeller):
    """
    MDSystem extends the Modeller class, providing useful methods for 
    preparing molecular systems for MD simulations in OpenMM.
    """
    def __init__(self, topology, positions, forcefield):
        Modeller.__init__(self, topology, positions)
        self.forcefield = forcefield
        self.simulation = None
        
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

    def save(self, pdbfile):
        """
        Save a PDB file representing this molecular system.

        Parameters
        ----------
        pdbfile : str
            PDB filename to which molecular system will be written
        """
        with open(pdbfile, "w") as outfile:
            PDBFile.writeFile(self.topology, self.positions, outfile)
        return
        
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
        
    def buildSimulation(self, integrator=LangevinIntegrator, dt=0.002*picoseconds,
                        temperature=298.15*kelvin, ensemble="NPT", posre=False,
                        posre_sel="not water and not (element Na or element Cl) and not element H",
                        efx=False, ef=(0,0,0), ef_sel="all", nonbondedMethod=PME,
                        nonbondedCutoff=1.*nanometer, constraints=HBonds, rigidWater=True, exceptions=[],
                        filePrefix="traj", saveTrajectory=False, trajInterval=500,
                        saveStateData=False, stateDataInterval=250):
        """
        Build a simulation context from the system. The simulation is
        then available as an attribute.
        """

        # If simulation exists, close any reporters
        if self.simulation is not None:
            for reporter in self.simulation.reporters:
                try:
                    reporter.close()
                except:
                    continue

        # Build system
        system = self.forcefield.createSystem(self.topology, nonbondedMethod=nonbondedMethod, 
                                              nonbondedCutoff=nonbondedCutoff, 
                                              constraints=constraints, rigidWater=rigidWater)

        # Setup MD simulation
        integrator = integrator(temperature, 1/picosecond, dt)

        # Add position restraints that can be tapered off during simulation
        if posre:
            force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
            force.addGlobalParameter("k", 5.0*kilocalories_per_mole/angstroms**2)
            force.addPerParticleParameter("x0")
            force.addPerParticleParameter("y0")
            force.addPerParticleParameter("z0")
            for i in self.select(posre_sel):
                force.addParticle(int(i), self.positions[i].value_in_unit(nanometers))
            system.addForce(force)

        # Add external electric field
        if efx:
            force = CustomExternalForce('(Ex*charge*x)+(Ey*charge*y)+(Ez*charge*z)')
            force.addGlobalParameter("Ex", ef[0])
            force.addGlobalParameter("Ey", ef[1])
            force.addGlobalParameter("Ez", ef[2])
            force.addPerParticleParameter("charge")
            es_forces = system.getForce(3)
            system.addForce(force)
            for i in self.select(ef_sel):
                i = int(i)
                charge = es_forces.getParticleParameters(i)[0]
                force.addParticle(i, [charge])

        # Setup exceptions in nonbonded forces if provided
        nonbonded = system.getForce(3)
        for atom1, atom2 in exceptions:
            nonbonded.addException(int(atom1), int(atom2), 0.0, 0.0, 0.0, True)
            
        # Setup barostat for NPT ensemble
        if ensemble == "NPT":
            barostat = MonteCarloBarostat(1.0*bar, temperature, 25)
            system.addForce(barostat)

        # Add simulation
        self.simulation = Simulation(self.topology, system, integrator)

        # Initialize particle positions
        self.simulation.context.setPositions(self.positions)

        # Add reporters
        if saveTrajectory:
            self.simulation.reporters.append(HDF5Reporter(f"{filePrefix}.h5", trajInterval))
        if saveStateData:
            self.simulation.reporters.append(StateDataReporter(f"{filePrefix}.csv", stateDataInterval, step=True, time=True, volume=True, totalEnergy=True, temperature=True, elapsedTime=True))
        
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
        self.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        
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

        # If simulation step is 0, output the starting configuration to
        # reporters
        if self.simulation.currentStep == 0:
            for reporter in self.simulation.reporters:
                report = reporter.describeNextReport(self.simulation)
                state  = self.simulation.context.getState(*report[1:])
                reporter.report(self.simulation, state)
                
        # Update positions
        state = self.simulation.context.getState(getPositions=True)
        self.positions = state.getPositions()
        self.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        
        return self

    def equilibrate(self, simtime=1.*nanoseconds, temperature=300*kelvin, posre=True):
        """
        Minimizes and equilibrate an MDSystem object. If position restraints
        are applied, it will taper the restraints over the course of the 
        simulation. This method assumes that MDSystem.buildSimulation() has
        already been called.

        Parameters
        ----------
        simtime : simtk.unit
            Total simulation time to use for equilibration
        temperature : OpenMM.unit.Quantity(unit=kelvin)
            Temperature to use to initialize velocities
        posre : bool
            If True, position restraints have been applied to simulation object
        """
        return equilibrate.equilibrate(self, simtime, temperature, posre)

    def calmdown(self, posre=True):
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
        posre : bool
            If true, position restraints are applied to non-water, 
            non-ion heavy atoms in the system.
        """
        return calmdown.calmdown( self, posre=posre)

    def getCharges(self, selection):
        """
        Get partial charges associated with atoms in selections. 

        Parameters
        ----------
        selection : str
            MDTraj-style atom selection string

        Returns
        -------
        np.ndarray
            Partial charges assigned by forcefield to atoms selected by
            selection string
        """
        remove = False
        if self.simulation is None:
            self.buildSimulation()
            remove = True

        force = self.simulation.system.getForce(3)
        indices = self.select(selection)
        charges = [ force.getParticleParameters(int(i))[0].value_in_unit(elementary_charge) for i in indices ]

        if remove:
            self.simulation = None
        
        return np.array(charges)
