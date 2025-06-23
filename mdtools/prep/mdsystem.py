"""
mdsystem.py: Extends Modeller class from OpenMM with useful methods

Author: Jack Greisman <greisman@g.harvard.edu>
        Ziyuan Zhao   <ziyuanzhao@fas.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.1"

from simtk.unit import *
from openmm import *
from openmm.app import *
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
        self._NonbondedForceIndex = None
        self._more_reporters = False

    def _getIndexOfNonbondedForce(self, system=None):
        if not system:
            system = self.simulation.system
        if not self._NonbondedForceIndex:
            for i, force in enumerate(system.getForces()):
                if isinstance(force, NonbondedForce):
                    self._NonbondedForceIndex = i
                    break
        return self._NonbondedForceIndex        

        
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
        
    def buildSimulation(self, integrator=LangevinMiddleIntegrator, dt=0.002*picoseconds,
                        temperature=298.15*kelvin, ensemble="NPT", posre=False,
                        posre_sel="not water and not (element Na or element Cl) and not element H",
                        efx=False, ef=(0,0,0), ef_sel="all", nonbondedMethod=PME,
                        nonbondedCutoff=1.*nanometer, constraints=HBonds, rigidWater=True, exceptions=[],
                        filePrefix="traj", saveTrajectory=False, trajInterval=500, saveVelocities=False,
                        saveStateData=False, stateDataInterval=250, atomSubset=None, thermalize=True,
                        hydrogenMass=1*amu, reporters=None):
        """Build a simulation context from the system. The simulation is
        then available as an attribute.

        Parameters
        ----------
        integrator : Openmm Integrator, optional
            Integrator for computing the MD trajectory, by default LangevinMiddleIntegrator
        dt : simtk.Unit, optional
            Time step size for integration, by default 0.002*picoseconds
        temperature : simtk.Unit, optional
            Temperature for the thermostat, by default 298.15*kelvin
        ensemble : str, optional
            Statistical ensemble for the simulation, by default "NPT"
        posre : bool, optional
            Whether to apply position restraints (posre), by default False
        posre_sel : str, optional
            Rule for selecting atoms to apply the posre, by default 
            "not water and not (element Na or element Cl) and not element H"
        efx : bool, optional
            Whether to apply electric field (EF) during simulation, by default False
        ef : tuple, optional
            Direction of the uniform EF, by default (0,0,0)
        ef_sel : str, optional
            Rule for selecting atoms that will feel the EF, by default "all"
        nonbondedMethod : OpenMM NonbondedForce, optional
            Model for nonbonded interactions between particles, by default PME
        nonbondedCutoff : simtk.Unit, optional
            Cutoff distance for nonbonded interactions, by default 1.*nanometer
        constraints : OpenmM Constraints, optional
            Constraints used for the simulation, by default HBonds
        rigidWater : bool, optional
            Whether to treat water molecules as rigid (e.g., 3-pt models), by default True
        exceptions : list, optional
            List of atoms that will be excluded from nonbonded forces treatment, 
            by default []
        filePrefix : str, optional
            Prefix to the saved files, by default "traj"
        saveTrajectory : bool, optional
            Whether to save trajectory from simulation, by default False
        trajInterval : int, optional
            Frequency of saving trajectory specified as the number of frames in between, by default 500
        saveVelocities : bool, optional
            Whether to save velocity from simulation, by default False
        saveStateData : bool, optional
            Whether to save other state data from simulation, by default False
        stateDataInterval : int, optional
            Frequency of saving other state data, by default 250
        atomSubset : Any | None, optional
            Indices of the subset of atoms to record trajectory for, 
            if None, all atoms will be recorded, by default None
        thermalize : bool, optional
            If True, initialize velocities according to Maxwell-Boltzmann 
            distribution, by default True
        hydrogenMass : simtk.Unit, optional
            Hydrogen mass, by default 1*amu
        reporters : List, optional
            List of reporters for collecting any additional information, by default None.
            A reporter is defined as a 5-tuple of (filePrefix, offset, trajInterval, 
            stateDataInterval, atomSubset)

        Returns
        -------
        mdtools.MDSystem
            Returns self, a modified MD system 
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
                                              constraints=constraints, rigidWater=rigidWater,
                                              hydrogenMass=hydrogenMass)

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

        # Add external electric field (specified as potential energy)
        if efx:
            force = CustomExternalForce('(-1*Ex*charge*x)+(-1*Ey*charge*y)+(-1*Ez*charge*z)')
            force.addGlobalParameter("Ex", ef[0])
            force.addGlobalParameter("Ey", ef[1])
            force.addGlobalParameter("Ez", ef[2])
            force.addPerParticleParameter("charge")
            es_forces = system.getForce(self._getIndexOfNonbondedForce(system))
            system.addForce(force)
            for i in self.select(ef_sel):
                i = int(i)
                charge = es_forces.getParticleParameters(i)[0]
                force.addParticle(i, [charge])

        # Setup exceptions in nonbonded forces if provided
        nonbonded = system.getForce(self._getIndexOfNonbondedForce(system))
        for atom1, atom2 in exceptions:
            nonbonded.addException(int(atom1), int(atom2), 0.0, 0.0, 0.0, True)
            
        # Setup barostat for NPT ensemble
        if ensemble == "NPT":
            barostat = MonteCarloBarostat(1.0*bar, temperature, 25)
            system.addForce(barostat)

        # Add simulation
        self.simulation = Simulation(self.topology, system, integrator)
        
        # Initialize particle positions and velocities
        self.simulation.context.setPositions(self.positions)
        if thermalize:
            self.thermalize(temperature=temperature)
        
        # Add reporters
        # We extend its functionality to allow adding multiple reporters,
        # and each reporter can have a different offset when it starts collecting data
        if reporters is None:
            if saveTrajectory:
                self.simulation.reporters.append(HDF5Reporter(f"{filePrefix}.h5", trajInterval, atomSubset=atomSubset, velocities=saveVelocities))
            if saveStateData:
                self.simulation.reporters.append(StateDataReporter(f"{filePrefix}.csv", stateDataInterval, step=True, time=True, volume=True, totalEnergy=True, temperature=True, elapsedTime=True))
            self._more_reporters = False
        else: # we assume it is a list of tuples (filePrefix, offset, trajInterval, stateDataInterval, atomSubset)
            self.reporters = reporters
            self._more_reporters = True

        return self

    def thermalize(self, temperature, randomSeed=None):
        """
        Set velocities of all particles to random values chosen from a 
        Maxwell-Boltzmann distribution for the given temeprature. 

        Parameters
        ----------
        temperature : float
            Temperature for which to sample velocities (Kelvin)
        randomSeed :  int
            Seed for random number generator
        """
        if randomSeed:
            return self.simulation.context.setVelocitiesToTemperature(temperature, randomSeed)
        else:
            return self.simulation.context.setVelocitiesToTemperature(temperature)

    def saveCheckpoint(self, filename):
        """
        Save a checkpoint of the simulation to a file.

        Parameters
        ----------
        filename : str
            File to which checkpoint will be saved
        """
        return self.simulation.saveCheckpoint(filename)

    def loadCheckpoint(self, filename):
        """
        Load a checkpoint of the simulation from a file.

        Parameters
        ----------
        filename : str
            File from which checkpoint will be loaded
        """
        return self.simulation.loadCheckpoint(filename)
    
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
        if isinstance(time, int):
            return time
        else:
            chemtime = time.in_units_of(picoseconds)
            dt = self.simulation.integrator.getStepSize()
            return int(np.ceil(chemtime / dt))
        
    def simulate(self, n, outputStartingFrame=True, reportLargeForceThreshold=-1):
        """
        Simulate the system for the given number of steps. If n is a 
        simtk.Unit of time, the number of steps are chosen to simulate
        for the indicated chemical time. 

        Parameters
        ----------
        n : int or simtk.unit
            Number of steps or chemical time of simulation
        outputStartingFrame : bool
            Whether to output the initial frame of a simulation
        reportLargeForceThreshold : int
            If <= 0, will not report; otherwise, print a list of
            all atoms with net forces on them exceeding the
            threshold in magnitude
        """

        # If simulation step is 0, output the starting configuration
        if self.simulation.currentStep == 0 and outputStartingFrame:
            for reporter in self.simulation.reporters:
                report = reporter.describeNextReport(self.simulation)
                state  = self.simulation.context.getState(*report[1:])
                reporter.report(self.simulation, state)
        n = self._time2steps(n)
        while n > 0:
            if self._more_reporters:
                next_offset = self._time2steps(self.reporters[0][1])
                n_steps = min(n, next_offset - self.simulation.currentStep)
                self.simulation.step(n_steps)
                print(next_offset, self.simulation.currentStep, n_steps)
                # Append new reporter to the list of reporters attached to the simulation if necessary
                filePrefix, offset, trajInterval, stateDataInterval, atomSubset = self.reporters.pop(0)
                # Requires rebuilding mdtraj
                if trajInterval > 0:
                    self.simulation.reporters.append(HDF5Reporter(f"{filePrefix}.h5", trajInterval, atomSubset=atomSubset, startTime=self.simulation.currentStep))
                # This doesn't work as it currently stands!
                if stateDataInterval > 0:
                    self.simulation.reporters.append(StateDataReporter(f"{filePrefix}.csv", stateDataInterval, step=True, time=True, volume=True, totalEnergy=True, temperature=True, elapsedTime=True))
                if len(self.reporters) == 0:
                    self._more_reporters = False
                n -= n_steps
            else:
                self.simulation.step(n)
                n = 0

        # Optionally report large forces
        if reportLargeForceThreshold > 0:
            state = self.simulation.context.getState(getForces=True)
            netforces = np.linalg.norm(state.getForces(asNumpy=True), axis=1)
            indices = np.where(np.isnan(netforces) | (netforces > reportLargeForceThreshold))[0]
            atoms = list(self.topology.atoms())
            print("The following atoms experience large net forces exceeding the threshold", reportLargeForceThreshold)
            [print(f"{atoms[idx]}, net F = {netforces[idx]}") for idx in indices]

        # Update positions
        state = self.simulation.context.getState(getPositions=True)
        self.positions = state.getPositions()
        self.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        
        return self

    def equilibrate(self, simtime=1.*nanoseconds, temperature=300*kelvin, posre=True, reportLargeForceThreshold=-1):
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
        reportLargeForceThreshold : int
            If <= 0, will not report; otherwise, print a list of
            all atoms with net forces on them exceeding the
            threshold in magnitude
        """
        return equilibrate.equilibrate(self, simtime, temperature, posre, reportLargeForceThreshold=reportLargeForceThreshold)

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

        force = self.simulation.system.getForce(self._getIndexOfNonbondedForce())
        indices = self.select(selection)
        charges = [ force.getParticleParameters(int(i))[0].value_in_unit(elementary_charge) for i in indices ]

        if remove:
            self.simulation = None
        
        return np.array(charges)
