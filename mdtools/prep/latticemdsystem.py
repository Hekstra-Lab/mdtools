"""
latticemdsystem.py: Provides LatticeMDSystem class for preparing MD
                    simulations of crystal lattices

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import mdtraj
from .mdsystem import MDSystem
import numpy as np
from math import floor
from scipy.spatial.distance import cdist
from itertools import product
from .xray import sg_canonicalize, sg_sym_to_mat_list
from mdtools.equilibration import squeeze

class LatticeMDSystem(MDSystem):
    """
    LatticeMDSystem provides methods for preparing MD simulations of 
    crystal lattices in OpenMM.
    """
    def __init__(self, topology, positions, forcefield, spacegroup):
        MDSystem.__init__(self, topology, positions, forcefield)
        self.spacegroup = spacegroup

    def addSolvent(self, distThreshold=3.2*angstroms, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True):
        """
        Solvate the unit cell. Build a waterbox that does not clash with the
        existing solute 
        """
        # Build up waterbox
        vectors = self.topology.getPeriodicBoxVectors()
        waterbox = self._initializeLatticeMDSystem()
        super(LatticeMDSystem, waterbox).addSolvent(self.forcefield, boxVectors=vectors, neutralize=False)
        waterbox._imageUnitCell(preserveMolecules=True)
        waterOs = [ a.index for a in waterbox.topology.atoms()
                    if a.element.symbol != "H" ]
        waterOs = np.array(waterOs)
        waterpos = np.array(waterbox.positions.value_in_unit(nanometers))
        waterpos = waterpos[waterOs]

        # Precompute heavy atom positions in 27 neighboring periodic images
        dims = internal.unitcell.computeLengthsAndAngles(vectors)
        heavyatoms = [ a.index for a in self.topology.atoms()
                       if a.element.symbol != "H" ]
        pos = np.array(self.positions.value_in_unit(nanometers))
        pos = pos[heavyatoms]
        pos = pos[np.newaxis, :, :]
        neighborhood = np.repeat(pos, 27, axis=0)
        cells = [-1, 0, 1]
        transmat = np.array([ uc for uc in product(cells, cells, cells) ])*dims[:3]
        neighborhood += transmat.reshape(27, 1, 3)

        # Compute all distances between water and unitcell atoms
        dists = [ cdist(waterpos, neighbor).min(axis=1) for neighbor in neighborhood ]

        # Find and remove clashing waters
        shortest = np.array(dists).min(axis=0)
        clashes = np.where(shortest < distThreshold.value_in_unit(nanometers))[0]
        clashers = waterOs[clashes]
        clashingWaters = [ a.residue for a in waterbox.topology.atoms() if a.index in clashers]
        waterbox.delete(clashingWaters)

        # Neutralize system
        if neutralize:

            # Determine whether divalent cation
            divalentCation = "2" in positiveIon
            
            # Determine ion type
            posIonElements = {'Cs+':element.cesium,
                              'K+':element.potassium,
                              'Li+':element.lithium,
                              'Na+':element.sodium,
                              'Rb+':element.rubidium,
                              'Mn2+':element.manganese,
                              'Mg2+':element.magnesium,
                              'Ca2+':element.calcium}
            negIonElements = {'Cl-':element.chlorine,
                              'Br-':element.bromine,
                              'F-':element.fluorine,
                              'I-':element.iodine}
            if positiveIon not in posIonElements:
                raise ValueError('Illegal value for positive ion: %s' % positiveIon)
            if negativeIon not in negIonElements:
                raise ValueError('Illegal value for negative ion: %s' % negativeIon)
            positiveElement = posIonElements[positiveIon]
            negativeElement = negIonElements[negativeIon]
        
            # Determine how many ions to add
            system = self.forcefield.createSystem(self.topology)
            nonbonded = None
            for i in range(system.getNumForces()):
                if isinstance(system.getForce(i), NonbondedForce):
                    nonbonded = system.getForce(i)
            if nonbonded is None:
                raise ValueError("The ForceField does not specify a NonbondedForce")
            totalCharge = int(floor(0.5+sum([ nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge) for i in range(system.getNumParticles()) ])))
            numIons = waterbox.topology.getNumResidues()*ionicStrength/(55.4*molar) # Pure water is about 55.4 molar (depending on temperature)
            numPairs = int(floor(numIons+0.5))

            # Add ions by replacing random waters
            ions = Modeller(Topology(), Quantity([]*nanometer))
            ions.topology.setPeriodicBoxVectors(vectors)
            newChain = ions.topology.addChain("ions")

            while totalCharge != 0:
                if totalCharge < 0:
                    ions, waterbox = self._waterToIon(ions, waterbox, positiveElement)
                    if divalentCation:
                        totalCharge += 2
                    else:
                        totalCharge += 1
                else:
                    ions, waterbox = self._waterToIon(ions, waterbox, negativeElement)
                    totalCharge -= 1

            for i in range(numPairs):
                ions, waterbox = self._waterToIon(ions, waterbox, positiveElement)
                ions, waterbox = self._waterToIon(ions, waterbox, negativeElement)
                # If divalent cation, neutralize with 2 negative ions
                if divalentCation:
                    ions, waterbox = self._waterToIon(ions, waterbox, negativeElement)

        self.add(waterbox.topology, waterbox.positions)
        if neutralize:
            self.add(ions.topology, ions.positions)

        # Image the unit cell
        self._imageUnitCell(preserveMolecules=True)
        return
            
    def _waterToIon(self, ions, waterbox, element):

        # Pick a water at random to replace with an ion
        waters = [ w for w in waterbox.topology.residues() ]
        waterpos = np.array(waterbox.positions.value_in_unit(nanometers))
        water = np.random.choice(waters)
        waterO = [ a for a in water.atoms() if a.element.symbol == "O"][0]
        waterbox.delete([water])

        # Place an ion at the position of the water oxygen
        chain = next(ions.topology.chains())
        newResidue = ions.topology.addResidue(element.symbol.upper(), chain)
        ions.topology.addAtom(element.symbol, element, newResidue)
        newPositions = ions.positions.value_in_unit(nanometers)
        newPositions.append(waterpos[waterO.index])
        ions.positions = newPositions*nanometers

        return ions, waterbox

    def _imageUnitCell(self, preserveMolecules=True):
        """
        Move molecules to best fit in the 000 unit cell
        
        Algorithm:
        1) Find molecules based on connected components of topology
        2) Image molecule if the center of geometry of the molecule lies in
        a neighboring unit cell
        """
        # Compute the de-orthogonalization matrix and get unit cell dimensions
        vectors = self.topology.getPeriodicBoxVectors()
        dims = list(internal.unitcell.computeLengthsAndAngles(vectors))
        lengths = np.array(dims[:3])
        orthmat = np.matrix(vectors.value_in_unit(nanometer))
        deorthmat = orthmat.I
        
        # Get positions of atoms in unit cell
        pos = np.array(self.positions.value_in_unit(nanometer))

        if preserveMolecules:
        
            # Find molecules in unit cell
            mols = self.findMolecules()

            # Image each molecule and move to 000 unit cell
            for mol in mols:

                indices = [ atom.index for atom in mol ]
                pos[indices] += self._imageAtoms(pos[indices].mean(axis=0), deorthmat, lengths)

        else:

            pos += self._imageAtoms(pos, deorthmat, lengths)
    
        self.positions = pos*nanometer
        return
    
    def _imageAtoms(self, pos, deorthmat, unitCellLengths):
        """
        Compute the image of the given atom positions that falls in the 000
        unit cell. Returns the translation matrix that moves the 
        corresponding atoms correctly.
        """
        fpos = np.array(pos*deorthmat)
        trans = unitCellLengths * ((fpos % 1) - fpos)
        return trans

    def _getUnitCellBasis(self):
        """
        For the unit cell specified in the Modeller object, compute the
        basis transformation as a 4x4 matrix
        """
        basis = np.zeros((4, 4))
        uc_vectors = self.topology.getPeriodicBoxVectors()
        uc_vectors = np.array(uc_vectors.value_in_unit(angstrom))
        basis[:3, :3] = uc_vectors.T
        basis[-1, -1] = 1.0
        return np.matrix(basis)

    def _initializeLatticeMDSystem(self):
        """Initialize an empty LatticeMDSystem object"""
        top = Topology()
        pos = Quantity((), angstroms)
        return LatticeMDSystem(top, pos, self.forcefield, self.spacegroup)

    def buildSuperCell(self, a=1, b=1, c=1, inplace=True):
        """
        Build a supercell using crystallographic symmetry
        
        Parameters
        ----------
        a : int
            Number of unit cell repeats along the a-axis
        b : int
            Number of unit cell repeats along the b-axis
        c : int
            Number of unit cell repeats along the c-axis
        """
        # Initialize Modeller object for buliding supercell
        supercell = self._initializeLatticeMDSystem()

        # Set periodic box vectors for supercell
        vectors = self.topology.getPeriodicBoxVectors()
        uc_dimensions = list(internal.unitcell.computeLengthsAndAngles(vectors))
        uc_dimensions[0] *= a
        uc_dimensions[1] *= b
        uc_dimensions[2] *= c
        newvectors = internal.unitcell.computePeriodicBoxVectors(*uc_dimensions)
        supercell.topology.setPeriodicBoxVectors(newvectors)
        
        # Build each unit cell
        for uc in product(range(a), range(b), range(c)):
            unitcell = self.buildUnitCell(*uc, inplace=False)
            supercell.add(unitcell.topology, unitcell.positions)

        if inplace:
            self.topology = supercell.topology
            self.positions = supercell.positions
        else:
            return supercell

    def buildUnitCell(self, a=0, b=0, c=0, inplace=True):
        """
        Build a supercell using crystallographic symmetry
        
        Parameters
        ----------
        a : int
            Index of unit cell along the a-axis
        b : int
            Index of unit cell along the b-axis
        c : int
            Index of unit cell along the c-axis
        inplace : bool
            Determines whether the unit cell is built in place or returned
            as a new LatticeMDSystem
        """
        # Get symmetry information
        sg = sg_canonicalize(self.spacegroup)
        matrices = sg_sym_to_mat_list(sg)
        basis = self._getUnitCellBasis()
        
        # Determine unit cell center
        pos = np.array(self.positions.value_in_unit(angstrom))
        top = self.getTopology()
        center = (pos.max(axis=0) + pos.min(axis=0))* 0.5
        center = np.matrix(center.tolist() + [1.0]).T
        center_cell = basis.I * center
        
        # Determine shift based on unit cell indices
        extra_shift = [ [float(i)] for i in (a, b, c) ]

        # Make Modeller object to build up unit cell
        unitcell = self._initializeLatticeMDSystem()
        unitcell.topology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
    
        # Build unit cell from asymmetric unit
        for i, mat in enumerate(matrices):

            # Compute TTT matrix
            mat = np.matrix(mat)
            shift = np.floor(mat * center_cell)
            mat[:3, 3] -= shift[:3, 0]
            mat[:3, 3] += extra_shift
            mat = basis * mat * basis.I
            
            # Separate TTT matrix components
            pretransmat = np.array(mat[3, :3].T).reshape(-1)
            rotmat = np.array(mat[:3, :3])
            posttransmat = np.array(mat[:3, 3]).reshape(-1)
            
            # Compute new positions
            newpos = pos.copy()
            newpos += pretransmat
            newpos = np.dot(rotmat, newpos.T).T
            newpos += posttransmat
            unitcell.add(top, newpos*angstrom)

        if inplace:
            self.topology = unitcell.topology
            self.positions = unitcell.positions
        else:
            return unitcell

    def squeeze(self, tolerance=0.003, maxIterations=10, maxSimtime=10*nanoseconds):
        return squeeze.squeeze(self, tolerance=tolerance,
                               maxIterations=maxIterations,
                               maxSimtime=maxSimtime)
