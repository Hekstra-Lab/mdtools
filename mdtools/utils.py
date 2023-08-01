from simtk.unit import *
import mdtraj
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv, norm, det
import gemmi
from typing import List
import reciprocalspaceship as rs
import sys, os, getopt, subprocess
###
from mdtraj.formats.hdf5 import *
from mdtraj.utils import in_units_of
from mdtraj.core.trajectory import Trajectory
import itertools

def save_hdf5(traj, filename, force_overwrite=True, mode='w'):
    """Save trajectory in HDF5 format. Retains crystallographic cell length and angles.

    Parameters
    ----------
    traj: mdtraj.Trajectory
        Input trajectory to save
    filename : _type_
        Name of the file for saving the trajectory
    force_overwrite : bool, optional
        If true, will overwrite existing content in the file, by default True
    mode : str, optional
        File open mode, by default 'w'
    """
    with HDF5TrajectoryFile(filename, mode, force_overwrite=force_overwrite) as f:
        f.write(coordinates=in_units_of(traj.xyz, Trajectory._distance_unit, f.distance_unit),
                time=traj.time,
                cell_lengths=in_units_of(traj.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                cell_angles=traj.unitcell_angles)
        f.topology = traj.topology

def _save_traj(traj, fname):
    """
    Wrapper function over `save_hdf5` to be used in this utils file
    """
    save_hdf5(traj, fname, force_overwrite=False, mode='a')

def getFieldStrength(e):
    """
    Convert from given unit of electric field strength to a value
    in OpenMM's standard unit
    """
    def convert(v):
        return (v*AVOGADRO_CONSTANT_NA).value_in_unit(kilojoule_per_mole / elementary_charge / nanometer)

    if isinstance(e, list):
        return [convert(e_c) for e_c in e]
    else:
        return convert(e)


def compute_RMSD_by_chain(target_traj, top, rule, chain_id_list = None):
    """Compute RMSD for each chain in a molecular system.

    Parameters
    ----------
    target_traj : mdtraj.Trajectory
        Trajectory of the molecular system of interest
    top : mdtraj.Topology
        Topology of the system
    rule : str
        mdtraj selection rule for chosing subset of atoms in each chain
    chain_id_list : List[int], optional
        List of chain_ids so RMSD values are only computed for these chains. 
        If None, will compute for all chains, by default None

    Returns
    -------
    Tuple (np.ndarray, np.ndarray, np.ndarray)
        Returns a 2D array of RMSD for each chain over time as the first item.
        The second and third items are the mean and std of RMSD across all chains.
    """
    chains_rmsd = []
    if chain_id_list == None:
        chain_id_list = range(top.n_chains)
    for chain_id in chain_id_list:
        atoms_selection = top.select("chainid " + str(chain_id) + " and " + rule)
        if(len(atoms_selection) == 0): # no valid atoms bound to a.a. residues
            continue
        chain_rmsd = mdtraj.rmsd(target_traj, target_traj, 0, atoms_selection)
        chains_rmsd.append(chain_rmsd)

    chains_rmsd = np.array(chains_rmsd)
    return chains_rmsd, np.mean(chains_rmsd, axis = 0), np.std(chains_rmsd, axis = 0)

def plot_RMSD_by_chain(chains_rmsd: np.ndarray, 
                       rule: str, 
                       n_frames: int, 
                       dn_frames: int = 100):
    """Plot the RMSD computed from `compute_RMSD_by_chain`

    Parameters
    ----------
    chains_rmsd : np.ndarray
        2D array of RMSD for each chain over time
    rule : str
        mdtraj selection rule for chosing subset of atoms in each chain
    n_frames : int
        Number of frames to include in the plot. 
    dn_frames: int
        Number of frames between every two points on the plot, by default 100.
    """
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, dn_frames), np.arange(0, n_frames, dn_frames) / 10)
    plt.ylabel("RMSD (nm)")
    for chid, chain_rmsd in enumerate(chains_rmsd):
        plt.plot(chain_rmsd, label = "chain " + str(chid))

    plt.legend(loc = 'upper right')


def plot_RMSD_by_chain_stat(avg_rmsd: np.ndarray, 
                            std_rmsd: np.ndarray, 
                            rule: str, 
                            n_frames: int,
                            dn_frames: int = 100):
    """Plot the RMSD statistics (avg and std) computed from `compute_RMSD_by_chain`

    Parameters
    ----------
    avg_rmsd : np.ndarray
        1D array of average RMSD across chains over time
    std_rmsd : np.ndarray
        1D array of standard deviation of RMSD across chains over time
    rule : str
        mdtraj selection rule for chosing subset of atoms in each chain
    n_frames : int
        Number of frames to include in the plot.
    dn_frames: int
        Number of frames between every two points on the plot, by default 100.
    """

    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, dn_frames), 
               np.arange(0, n_frames, dn_frames) / 10)
    plt.ylabel("RMSD (nm)")
    plt.plot(avg_rmsd, label = "averaged")
    plt.plot(avg_rmsd + std_rmsd, label = "+1 sd")
    plt.plot(avg_rmsd - std_rmsd, label = "-1 sd")
    plt.legend(loc = 'upper right')

def unwrap_atom_axis(traj: mdtraj.Trajectory):
    """Unwrap Cartesian coordinates of the molecular system along atomic indices.

    In a periodic system, during molecular simulation, the coordinates are only 
    meaningful upto some modulo operation. Thus, molecule may cut through the 
    boundary and be wrapped back at the opposite face of the system, affecting 
    subsequent analyses. This function removes the wrapping artifact by making 
    sure the neighboring atoms are in proximity.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Input trajectory, will be modified in place.
    """
    L_arr = traj.unitcell_lengths[0]
    for i in range(3):
        L=L_arr[i]
        eps=L/2 #nm
        x0i = traj.xyz[:,0:1,i]
        dxi = traj.xyz[:,1:,i]-traj.xyz[:,:-1,i]
        dxi_corr=dxi-L*(np.sign(dxi)*((np.abs(dxi)-eps)>0))
        xi_corr = np.cumsum(np.concatenate((x0i, dxi_corr), axis=1), axis=1)
        traj.xyz[:,:,i] = xi_corr

def unwrap_time_axis(traj, selection='all'):
    """Unwrap Cartesian coordinates of molecular system along the time axis.

    In a periodic system, during molecular simulation, the coordinates are only 
    meaningful upto some modulo operation. Over time, the center of mass of a 
    molecule may drift out of one boundary so the molecule suddenly appears at
    the opposite boundary in subsequent frames. Such discontinuous jumps in 
    molecular coordinates also hampers analyses and this function removes the 
    artifact by making sure the coordinates from neighboring frames are close.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Input trajectory, will be modified in place.
    selection : List[int] | str, optional
        Atom selection for fixing coordinates, by default 'all'
    """
    L_arr = traj.unitcell_lengths[0]
    if selection == 'all':
        selection = np.arange(traj.xyz.shape[1])
    for i in range(3):
        L=L_arr[i]
        eps=L/2 #nm
        dxi = traj.xyz[1:,selection,i]-traj.xyz[:-1,selection,i]
        dxi_corr=dxi-L*(np.sign(dxi)*((np.abs(dxi)-eps)>0))
        xi_corr = np.cumsum(np.concatenate((traj.xyz[0:1,selection,i], dxi_corr)), axis=0)
        traj.xyz[:,selection,i] = xi_corr

DEN = 24 # why this? See https://gemmi.readthedocs.io/en/latest/symmetry.html for an explanation...

def real_space_transform(xyz: np.ndarray, 
                         op: gemmi.Op, 
                         unitcell_vecs: np.ndarray, 
                         wrap: bool = True,
                         fractional: bool = False) -> np.ndarray:
    """Transforms an array of real coordinates by a symmetry operator and returns the new array.

    Parameters
    ----------
    xyz : np.ndarray
        n * 3 array that stores the real space coordinates to be transformed
    op : gemmi.Op
        The symmetry operation to be applied
    unitcell_vecs : np.ndarray
        3 * 3 array with rows as unit cell vectors
    wrap : bool, optional
        If true, fractional coordinates will be wrapped to be inside the unit box, by default True
    fractional : bool, optional
        If True, returns fractional coordinate instead, by default False

    Returns
    -------
    np.ndarray
        n * 3 array containing real coordinates as transformed by op.
    """
    tran = np.expand_dims(np.array(op.tran) / DEN, axis=0)
    rot = np.expand_dims(inv(np.expand_dims(np.array(op.rot), axis=0)[0] / DEN), axis=0)
    deorthmat = np.array(unitcell_vecs)
    orthmat = np.expand_dims(inv(deorthmat), axis=0)
    if fractional is True:
        deorthmat = np.expand_dims(np.identity(3), axis=0)
    else:
        deorthmat = np.expand_dims(deorthmat, axis=0)
    if wrap is True:
        return ((xyz @ orthmat) @ rot + tran) % 1 @ deorthmat
    else:
        return ((xyz @ orthmat) @ rot + tran) @ deorthmat

def trajectory_transform(traj: mdtraj.Trajectory, 
                         op: gemmi.Op, inplace=True) -> None:
    """Performs a transformation applied to all atomic coordinates in each frame
    of the input trajectory according to symmetry operator `op`.

    Assumes constant unit cell vector for speed.
    Assumes CoM does't jump around near the boundary of the cell.
    (use unwrapping functions on the trajectory first if needed)

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Input trajectory that the transformation will be applied to, frame by frame.
    op : gemmi.Op
        Symmetry operator to be applied, contains translation and rotation parts.
    inplace : bool, optional
        If true, will transform trajectory inplace, else return a modified copy, 
        by default True

    Returns
    -------
    mdtraj.Trajectory | None
        New trajectory if `inplace` is True.
    """

    tran = np.expand_dims(np.array(op.tran) / DEN, axis=0)
    rot = np.expand_dims(inv(np.expand_dims(np.array(op.rot), axis=0)[0] / DEN), axis=0)
    deorthmat = traj.unitcell_vectors[0]
    orthmat = np.expand_dims(inv(deorthmat), axis=0)
    deorthmat = np.expand_dims(deorthmat, axis=0)
    com = mdtraj.compute_center_of_mass(traj[0])[0:1]
    frac_com = ((com @ orthmat) @ rot + tran)
    com_shift = -(frac_com // 1)
    new_xyz = ((traj.xyz @ orthmat) @ rot + tran + com_shift) @ deorthmat
    if inplace:
        traj.xyz = new_xyz
    else:
        return new_xyz

"""
    Performs a transformation applied to all atomic coordinates in each frame of the input trajectory to
    restore each chain back to the first chain (ASU), returns the new trajectory.
    Assumes that each chain doesn't move too much in the unit cell so the same symop can be applied across all frames.
    Assumes reference chain is the one with id = 0.
    :param traj: input trajectory that the transformation will be applied to.
    :param sg: space group code.
    :param chain_ids: a list of chains to perform the operation on.
    :return: None
    """

def trajectory_revert_to_asu(traj: mdtraj.Trajectory, 
                             sg: int, 
                             chain_ids: List[int], 
                             threshold=0.05):
    """Performs a transformation applied to all atomic coordinates in each frame
    of the input trajectory to restore each chain back to the first chain (ASU), 
    returns the new trajectory.

    Assumes that each chain doesn't move too much in the unit cell so the same 
    symop can be applied across all frames.
    Assumes reference chain is the one with id = 0.
    
    Parameters
    ----------
    traj : mdtraj.Trajectory
        Input trajectory that the transformation will be applied to.
    sg : int
        Space group code describing the symmetry of the system.
    chain_ids : List[int]
        A list of chains to perform the operation on.
    threshold : float, optional
        Decimal fraction threshold that determines how much CoM drift compared to 
        the cell dimension is allowed between the reference and symop reverted ASU, 
        by default 0.05.

    Returns
    -------
    mdtraj.Trajectory
        New trajectory with all chains reverted to the first ASU position.
    """

    frame0 = traj[0]
    top = traj.topology
    ucv = traj.unitcell_vectors[0]
    char_length = det(ucv) ** (1 / 3)
    ops = list(gemmi.find_spacegroup_by_number(sg).operations().sym_ops)
    ops_inv = [op.inverse() for op in ops]
    ref_sel = top.select(f"chainid 0")
    ref_com = mdtraj.compute_center_of_mass(frame0.atom_slice(ref_sel))[0]
    new_traj = traj.atom_slice(ref_sel)
    for id in chain_ids:
        if id == 0:  # reference
            continue
        top_sel = top.select(f"chainid {id}")
        com = mdtraj.compute_center_of_mass(frame0.atom_slice(top_sel))[0]  # (3,)
        subtraj = traj.atom_slice(top_sel)
        for op, op_inv in zip(ops, ops_inv):
            new_com = real_space_transform(com, op_inv, ucv)
            ratio = norm(new_com - ref_com) / char_length
            if ratio < threshold:
                trajectory_transform(subtraj, op_inv)
                print(f"Chain {id} corresponds to symop {op.triplet()}")
                break
        new_traj = new_traj.stack(subtraj, keep_resSeq=False)
    return new_traj

def align_and_split_by_chain(traj: mdtraj.Trajectory, 
                             output_name: str, 
                             unitcell_ref: mdtraj.Trajectory, 
                             asu_ref: mdtraj.Trajectory, 
                             chainwise_alignment=True,
                             asu_reversion: bool = True, 
                             sg: int = None, 
                             atom_selection: List[int] | None = None):
    """Align and split the peptide chains in the trajectory.

    We provide three modes of alignment that are useful in different scenarios:
    1. Alignment of entire cell content: useful for studying structural deformation
    of entire crystal simulation system as it retains the relative displacements of
    chains from their equilibrium position. This can be done with 
    `chainwise_alignment = False, asu_reversion = False`.
    2. Alignment by symops: this is option 1 plus an additional reversion of chains
    back to first ASU position via all symops of the space group. It facilitates
    studying the anisotropic drifts among the chains due to external pertubation
    during a crystal simulation. Uses parameter `chainwise_alignment = False, 
    asu_reversion = True`.
    3. Alignment by RMSD minimization: useful for studying internal motions in
    each chain because the relative displacement among chains are eliminated via
    superposing each chain to a reference ASU. Uses parameter 
    `chainwise_alignment = True, asu_reversion = True/False`. 

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory to perform alignment on and extract chains from
    output_name : str
        Output file name. File will be saved in the **current** directory
    unitcell_ref : mdtraj.Trajectory
        Unit cell reference structure, used if `chainwise_alignment = False`
    asu_ref : mdtraj.Trajectory
        Asymmetric unit reference structure, used if `chainwise_alignment = True`
    chainwise_alignment : bool, optional
        Switches between alignment modes, see above for explanation, by default True
    asu_reversion : bool, optional
        Switches between alignment modes, see above for explanation, by default True
    sg : int, optional
        Space group number, by default None
    atom_selection : List[int] | None, optional
        List of atom indices to use for alignment, if None, use all atoms. 
        By default None
    """
    traj = traj.remove_solvent()
    n_frames = traj.xyz.shape[0]
    top = traj.topology
    asu_ref = asu_ref[0]

    if chainwise_alignment:
        traj.unitcell_lengths = np.repeat(asu_ref.unitcell_lengths, axis=0, repeats=n_frames)
        traj.unitcell_angles = np.repeat(asu_ref.unitcell_angles, axis=0, repeats=n_frames)
        traj.unitcell_vectors = np.repeat(asu_ref.unitcell_vectors, axis=0, repeats=n_frames)
    else:
        traj.superpose(unitcell_ref)
        if asu_reversion:
            traj = trajectory_revert_to_asu(traj, sg=sg, chain_ids=np.arange(traj.n_chains))

    for chain_id in np.arange(traj.n_chains):
        subtraj = traj.atom_slice(top.select(f"chainid == {chain_id}"))
        if chainwise_alignment:
            subtraj.superpose(asu_ref, atom_indices=atom_selection)
        _save_traj(subtraj, f"{output_name}_subtraj_{chain_id}.h5")

def save_snapshots_from_traj(target_traj: mdtraj.Trajectory, 
                             output_name: str, 
                             frame_offset: int, 
                             d_frame: int):
    """Handy function to save snapshots from trajectory

    Parameters
    ----------
    target_traj : mdtraj.Trajectory
        Trajectory to save from 
    output_name : str
        Output file name. File will be saved to the **current** directory
    frame_offset : int
        Initial frame to slice from
    d_frame : int
        Number of frames between every two saved snapshots
    """

    for i, frame in enumerate(target_traj[frame_offset::d_frame]):
#         frame.superpose(ref_traj, atom_indices = top.select("is_backbone"),
#                         ref_atom_indices = top.select("is_backbone")).save_pdb(f"{output_name}_{i}.pdb")
        frame.save_pdb(f"{output_name}_{i}.pdb")

def batch_annotate_spacegroup(input_name: str, 
                              max_frame: int, 
                              sg: str):
    """Annotate spacegroup for a batch of PDB files

    Parameters
    ----------
    input_name : str
        Input file to read from. Will read file from **current** directory
    max_frame : int
        Number of frames to read
    sg : str
        Space group number
    """
    for frame_id in range(max_frame):
        struct = gemmi.read_pdb(f"{input_name}_{frame_id}.pdb")
        struct.spacegroup_hm = sg
        struct.write_pdb(f"{input_name}_{frame_id}.pdb")

def batch_fmodel(input_name: str, 
                 max_frame: int, 
                 resolution: int = 1.5,
                 phenix_command: str = 'source /usr/local/phenix-1.20.1-4487/phenix_env.sh; phenix.fmodel', 
                 output_path: str | None = None):
    """Compute electron density map for a batch of PDB files. 

    Note: must make sure the files have spacegroup annotations first.

    Parameters
    ----------
    input_name : str
        Input file name
    max_frame : int
        Number of frames to process
    resolution : int, optional
        Electron density map resolution, by default 1.5 (angstrom)
    phenix_command : str, optional
        Command for running phenix.fmodel from your machine, 
        by default 'source /usr/local/phenix-1.20.1-4487/phenix_env.sh; phenix.fmodel'
    output_path : str | None, optional
        Output file path, if None outputs to current path, by default None
    """
    for frame_id in range(max_frame):
            subprocess.run(('' if output_path is None else f'cd {output_path};') +\
                            phenix_command + f' {input_name}_{frame_id}.pdb high_resolution={resolution}',
                             shell=True, stdout=subprocess.DEVNULL)

def average_structure_factors(input_name: str, 
                              max_frame: int):
    """Average structure factors for a batch of .mtz files

    Parameters
    ----------
    input_name : str
        Input file name
    max_frame : int
        Number of frames to process
    """
    dataset = rs.read_mtz(f"{input_name}_0.pdb.mtz")
    n_reflections = dataset.shape[0]

    complex_reflections = np.zeros(n_reflections, dtype='complex128')

    for frame in range(max_frame):
        dataset = rs.read_mtz(f"{input_name}_{frame}.pdb.mtz")
        complex_reflections = complex_reflections * (1 - 1/(frame + 1)) + np.array([amp*np.exp(np.pi*phase/180 * 1j) for [amp, phase] in dataset.to_numpy()]) / (frame + 1)


    dataset[:] = np.stack([np.abs(complex_reflections), np.angle(complex_reflections) / np.pi * 180]).T
    dataset.infer_mtz_dtypes(inplace = True)
    dataset.write_mtz(f"{input_name}_avg.mtz")

def _generate_flexible_fun(f, reader, writer):
    """Transforms f(input)                                      -> retval
               to f_flex(generalized_input, generalized_output) -> None / retval'
        
        An experimental feature that allows wrapping a simple function with single
        input and single output to handle multiple inputs or outputs.

        (this is not very relevant to the purpose of our package though...)
    """
    def f_flex(input, output):
        if isinstance(input, list):
            if isinstance(output, list) and len(output) == len(input):
                [writer(f(reader(i)), o) for i, o in zip(input, output)]
            elif output is None:
                return [f(reader(i)) for i in input]
            else:
                writer([f(reader(i)) for i in input], output)
        else:
            if output is None:
                return f(input)
            else:
                writer(f(reader(input)), output)
    return f_flex

def compute_net_dipole_moment(partial_charges: np.ndarray, 
                              input=None, 
                              output=None):
    """Compute net dipole moment from trajectory.

    Parameters
    ----------
    partial_charges : np.ndarray
        1D array containing atomic partial charges for all atoms in the system.
    input : Any, optional
        Input file name, or a list of file name, by default None
    output : Any, optional
        Output file name, or a list of file name, or None, by default None
    """
    def core(traj, args=None, kwargs=None):
        com = mdtraj.compute_center_of_mass(traj)
        xyz_minus_com = traj.xyz - com[:,None,:]
        net_dipole = np.sum(xyz_minus_com * partial_charges[None, :, None], axis=1)
        return net_dipole
    def writer(arr, output_name):
        np.save(output_name, arr)
    def reader(input):
        if isinstance(input, str):
            return mdtraj.load(input)
        return input
    return (_generate_flexible_fun(core, reader, writer))(input, output)

def mtz_to_cartesian_arr(mtz: rs.Dataset):
    """Convert structure factors from complex to Cartesian representation
    """
    return np.array([amp*np.exp(np.pi*phase/180 * 1j) for [amp, phase] in mtz.to_numpy()])

def cartesian_arr_to_polar(arr: np.array):
    """Convert structure factors from Cartesian to complex representation
    """
    return np.stack([np.abs(arr), np.angle(arr) / np.pi * 180]).T

def compute_all_difference_maps(file_name: str, 
                                n_chains: int = 4, 
                                phases: List[str] = ['pos', 'zero', 'neg']):
    """Compute internal and ordinary difference maps using all chains' maps.

    Parameters
    ----------
    file_name : str
        Input file name
    n_chains : int, optional
        Number of chains in the system, by default 4
    phases : List[str], optional
        List of phases of EF perturbation, by default ['pos', 'zero', 'neg']
    """
    # calculate internal diff map for all phases and pairs
    for phase in phases:
        for pair in itertools.combinations(range(n_chains), 2):
            chain_id_1, chain_id_2 = pair
            dataset1 = rs.read_mtz(f'{file_name}_chainwise_{phase}_subtraj_{chain_id_1}_avg.mtz')
            dataset2 = rs.read_mtz(f'{file_name}_chainwise_{phase}_subtraj_{chain_id_2}_avg.mtz')
            dataset1[:] = cartesian_arr_to_polar(mtz_to_cartesian_arr(dataset1) - mtz_to_cartesian_arr(dataset2))
            dataset1.infer_mtz_dtypes(inplace = True)
            dataset1.write_mtz(f'{file_name}_diff_{phase}_{chain_id_1}_{chain_id_2}.mtz')

    # calculate ordinary diff map for all phases and chains:
    for phase_pair in itertools.combinations(phases, 2):
        phase_1, phase_2 = phase_pair
        for chain_id in range(n_chains):
            dataset1 = rs.read_mtz(f'{file_name}_chainwise_{phase_1}_subtraj_{chain_id}_avg.mtz')
            dataset2 = rs.read_mtz(f'{file_name}_chainwise_{phase_2}_subtraj_{chain_id}_avg.mtz')
            dataset1[:] = cartesian_arr_to_polar(mtz_to_cartesian_arr(dataset1) - mtz_to_cartesian_arr(dataset2))
            dataset1.infer_mtz_dtypes(inplace = True)
            dataset1.write_mtz(f'{file_name}_diff_{chain_id}_{phase_1}_{phase_2}.mtz')
