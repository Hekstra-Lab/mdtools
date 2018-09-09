from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from math import floor
from scipy.spatial.distance import cdist
from itertools import product
import xray

def _findMolecules(top):
    """
    Use the list of bonds to divide up the Topology's atoms into molecules.

    This code was adapted from the MDTraj module
    """
    # Make a list of every other atom to which each atom is connected.
    num_atoms = top.getNumAtoms()
    atom_bonds = [[] for i in range(num_atoms)]
    for atom1, atom2 in top.bonds():
        atom_bonds[atom1.index].append(atom2.index)
        atom_bonds[atom2.index].append(atom1.index)
        
    # This is essentially a recursive algorithm, but it is reformulated as a loop to avoid
    # stack overflows.  It selects an atom, marks it as a new molecule, then recursively
    # marks every atom bonded to it as also being in that molecule.
    atom_molecule = [-1]*num_atoms
    num_molecules = 0
    for i in range(num_atoms):
        if atom_molecule[i] == -1:

            # Start a new molecule.
            atom_stack = [i]
            neighbor_stack = [0]
            molecule = num_molecules
            num_molecules += 1
                
            # Recursively tag all the bonded atoms.
            while len(atom_stack) > 0:
                atom = atom_stack[-1]
                atom_molecule[atom] = molecule
                while neighbor_stack[-1] < len(atom_bonds[atom]) and atom_molecule[atom_bonds[atom][neighbor_stack[-1]]] != -1:
                    neighbor_stack[-1] += 1
                if neighbor_stack[-1] < len(atom_bonds[atom]):
                    atom_stack.append(atom_bonds[atom][neighbor_stack[-1]])
                    neighbor_stack.append(0)
                else:
                    del atom_stack[-1]
                    del neighbor_stack[-1]

    # Build the final output.
    molecules = [set() for i in range(num_molecules)]
    for atom in top.atoms():
        molecules[atom_molecule[atom.index]].add(atom.index)
    return molecules

def _solvateUnitCell(unitcell, forcefield, distThreshold=3.2*angstroms, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True):
    """
    Solvate the unit cell. Build a waterbox that does not clash with the
    existing solute 
    """
    # Build up waterbox
    vectors = unitcell.topology.getPeriodicBoxVectors()
    waterbox = Modeller(Topology(), Quantity([]*nanometers))
    waterbox.addSolvent(forcefield, boxVectors=vectors, neutralize=False)
    _imageUnitCell(waterbox, preserveMolecules=True)
    waterOs = [ a.index for a in waterbox.topology.atoms()
                if a.element.symbol != "H" ]
    waterOs = np.array(waterOs)
    waterpos = np.array(waterbox.positions.value_in_unit(nanometers))
    waterpos = waterpos[waterOs]

    # Precompute heavy atom positions in 27 neighboring periodic images
    dims = internal.unitcell.computeLengthsAndAngles(vectors)
    heavyatoms = [ a.index for a in unitcell.topology.atoms()
                   if a.element.symbol != "H" ]
    pos = np.array(unitcell.positions.value_in_unit(nanometers))
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

        # Determine ion type
        posIonElements = {'Cs+':element.cesium,
                          'K+':element.potassium,
                          'Li+':element.lithium,
                          'Na+':element.sodium,
                          'Rb+':element.rubidium}
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
        system = forcefield.createSystem(unitcell.topology)
        nonbonded = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), NonbondedForce):
                nonbonded = system.getForce(i)
        if nonbonded is None:
            raise ValueError("The ForceField does not specify a NonbondedForce")
        totalCharge = int(floor(0.5+sum([ nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge) for i in range(system.getNumParticles()) ])))
        print("Total System Charge: {}".format(totalCharge))
        numIons = waterbox.topology.getNumResidues()*ionicStrength/(55.4*molar) # Pure water is about 55.4 molar (depending on temperature)
        numPairs = int(floor(numIons+0.5))
        print("Number of Ion Pairs: {}".format(numPairs))

        # Add ions by replacing random waters
        ions = Modeller(Topology(), Quantity([]*nanometer))
        ions.topology.setPeriodicBoxVectors(vectors)
        newChain = ions.topology.addChain("ions")

        for i in range(abs(totalCharge)):
            if totalCharge < 0:
                ions, waterbox = _waterToIon(ions, waterbox, positiveElement)
            else:
                ions, waterbox = _waterToIon(ions, waterbox, negativeElement)

        for i in range(numPairs):
            ions, waterbox = _waterToIon(ions, waterbox, positiveElement)
            ions, waterbox = _waterToIon(ions, waterbox, negativeElement)

    unitcell.add(waterbox.topology, waterbox.positions)
    if neutralize:
        unitcell.add(ions.topology, ions.positions)

    return unitcell

def _waterToIon(ions, waterbox, element):

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
    
def _imageUnitCell(unitcell, preserveMolecules=True):
    """
    Move molecules to best fit in the 000 unit cell

    Algorithm:
    1) Find molecules based on connected components of topology
    2) Image molecule if the center of geometry of the molecule lies in
       a neighboring unit cell
    """
    # Compute the de-orthogonalization matrix and get unit cell dimensions
    vectors = unitcell.topology.getPeriodicBoxVectors()
    dims = list(internal.unitcell.computeLengthsAndAngles(vectors))
    lengths = np.array(dims[:3])
    orthmat = np.matrix(vectors.value_in_unit(nanometer))
    deorthmat = orthmat.I

    # Get positions of atoms in unit cell
    pos = np.array(unitcell.positions.value_in_unit(nanometer))

    if preserveMolecules:
        
        # Find molecules in unit cell
        mols = _findMolecules(unitcell.topology)

        # Image each molecule and move to 000 unit cell
        for mol in mols:

            indices = list(mol)
            pos[indices] += _imageAtoms(pos[indices].mean(axis=0), deorthmat, lengths)

    else:

        pos += _imageAtoms(pos, deorthmat, lengths)
    
    unitcell.positions = pos*nanometer
    return unitcell

def _imageAtoms(pos, deorthmat, unitCellLengths):
    """
    Compute the image of the given atom positions that falls in the 000
    unit cell. Returns the translation matrix that moves the 
    corresponding atoms correctly.
    """
    fpos = np.array(pos*deorthmat)
    trans = unitCellLengths * ((fpos % 1) - fpos)
    return trans

def _getUnitCellBasis(modeller):
    """
    For the unit cell specified in the Modeller object, compute the
    basis transformation as a 4x4 matrix
    """
    basis = np.zeros((4, 4))
    uc_vectors = modeller.topology.getPeriodicBoxVectors()
    uc_vectors = np.array(uc_vectors.value_in_unit(angstrom))
    basis[:3, :3] = uc_vectors.T
    basis[-1, -1] = 1.0
    return np.matrix(basis)

def _initializeModeller():
    """Initialize an empty Modeller object"""
    return Modeller(Topology(), Quantity((), angstroms))

def buildSuperCell(modeller, spacegroup, a=1, b=1, c=1):
    """
    Build a supercell using crystallographic symmetry

    Parameters
    ----------
    modeller : openmm.Modeller
        OpenMM Modeller object that will be tiled to build the
        supercell
    spacegroup : str
        Crystallographic space group symmetry to use to construct
        supercell
    a : int
        Number of unit cell repeats along the a-axis
    b : int
        Number of unit cell repeats along the b-axis
    c : int
        Number of unit cell repeats along the c-axis
    """
    # Initialize Modeller object for buliding supercell
    supercell = _initializeModeller()

    # Set periodic box vectors for supercell
    vectors = modeller.topology.getPeriodicBoxVectors()
    uc_dimensions = list(internal.unitcell.computeLengthsAndAngles(vectors))
    uc_dimensions[0] *= a
    uc_dimensions[1] *= b
    uc_dimensions[2] *= c
    newvectors = internal.unitcell.computePeriodicBoxVectors(*uc_dimensions)
    supercell.topology.setPeriodicBoxVectors(newvectors)
    
    # Build each unit cell
    for uc in product(range(a), range(b), range(c)):
        unitcell = buildUnitCell(modeller, spacegroup, *uc)
        supercell.add(unitcell.topology, unitcell.positions)
    
    return supercell

def buildUnitCell(modeller, spacegroup, a=0, b=0, c=0):
    """
    Build a supercell using crystallographic symmetry

    Parameters
    ----------
    modeller : openmm.Modeller
        OpenMM Modeller object that will be tiled to build the
        supercell
    spacegroup : str
        Crystallographic space group symmetry to use to construct
        supercell
    a : int
        Index of unit cell along the a-axis
    b : int
        Index of unit cell along the b-axis
    c : int
        Index of unit cell along the c-axis
    """
    # Get symmetry information
    sg = xray.sg_canonicalize(spacegroup)
    matrices = xray.sg_sym_to_mat_list(sg)
    basis = _getUnitCellBasis(modeller)
    
    # Determine unit cell center
    pos = np.array(modeller.positions.value_in_unit(angstrom))
    top = modeller.getTopology()
    center = (pos.max(axis=0) + pos.min(axis=0))* 0.5
    center = np.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center
    
    # Determine shift based on unit cell indices
    extra_shift = [ [float(i)] for i in (a, b, c) ]

    # Make Modeller object to build up unit cell
    unitcell = _initializeModeller()
    unitcell.topology.setPeriodicBoxVectors(modeller.topology.getPeriodicBoxVectors())
    
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
        newpos = np.dot(newpos, rotmat)
        newpos += posttransmat
        unitcell.add(top, newpos*angstrom)

    return unitcell
