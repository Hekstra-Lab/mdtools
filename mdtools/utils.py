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

def save_hdf5(self, filename, force_overwrite=True, mode='w'):
    with HDF5TrajectoryFile(filename, mode, force_overwrite=force_overwrite) as f:
                f.write(coordinates=in_units_of(self.xyz, Trajectory._distance_unit, f.distance_unit),
                        time=self.time,
                        cell_lengths=in_units_of(self.unitcell_lengths, Trajectory._distance_unit, f.distance_unit),
                        cell_angles=self.unitcell_angles)
                f.topology = self.topology

def _save_traj(traj, fname):
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
    chains_rmsd = []
    n_frames = target_traj.n_frames
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

def plot_RMSD_by_chain(chains_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    for chid, chain_rmsd in enumerate(chains_rmsd):
        plt.plot(chain_rmsd, label = "chain " + str(chid))

    plt.legend(loc = 'upper right')


def plot_RMSD_by_chain_stat(avg_rmsd, std_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    plt.plot(avg_rmsd, label = "averaged")
    plt.plot(avg_rmsd + std_rmsd, label = "+1 sd")
    plt.plot(avg_rmsd - std_rmsd, label = "-1 sd")
    plt.legend(loc = 'upper right')

def unwrap_atom_axis(traj):
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
def real_space_transform(xyz: np.ndarray, op: gemmi.Op, unitcell_vecs: np.ndarray, wrap: bool = True,
                         fractional: bool = False) -> np.ndarray:
    """
    Transforms an array of real coordinates by a symmetry operator and returns the new array.
    :param xyz: n * 3 array that stores the real space coordinates to be transformed
    :param op: the symop to be applied
    :param unitcell_vecs: 3 * 3 array with rows as unit cell vectors
    :param wrap: if True, fractional coordinates will be wrapped to be inside the unit box
    :param fractional: if True, returns fractional coordinate instead
    :return: n * 3 array containing real coordinates as transformed by op.
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

def trajectory_transform(traj: mdtraj.Trajectory, op: gemmi.Op, inplace=True) -> None:
    """
    Performs a transformation applied to all atomic coordinates in each frame of the input trajectory according to symmetry operator `op`
    Assumes constant unit cell vector for speed.
    Assumes CoM does't jump around near the boundary of the cell (use smart wrap on the trajectory first if needed)
    :param traj: Input trajectory that the transformation will be applied to, frame by frame.
    :param op: Symmetry operator to be applied, contains translation and rotation parts.
    :return: None.
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


def trajectory_revert_to_asu(traj: mdtraj.Trajectory, sg: int, chain_ids: List[int], threshold=0.05):
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
        # print(f"Looking at chain {id}, com {com}, refcom {ref_com}")
        for op, op_inv in zip(ops, ops_inv):
            new_com = real_space_transform(com, op_inv, ucv)
            ratio = norm(new_com - ref_com) / char_length
            # print(f"Looking at symop {op.triplet()}, ratio {ratio}")
            if ratio < threshold:
                trajectory_transform(subtraj, op_inv)
                print(f"Chain {id} corresponds to symop {op.triplet()}")
                break
        new_traj = new_traj.stack(subtraj, keep_resSeq=False)
    return new_traj

def align_and_split_by_chain(traj, output_name, unitcell_ref, asu_ref, chainwise_alignment=True,
                                   asu_reversion=True, sg=None, atom_selection=None):
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

def save_snapshots_from_traj(target_traj, output_name, frame_offset, d_frame):
    for i, frame in enumerate(target_traj[frame_offset::d_frame]):
#         frame.superpose(ref_traj, atom_indices = top.select("is_backbone"),
#                         ref_atom_indices = top.select("is_backbone")).save_pdb(f"{output_name}_{i}.pdb")
        frame.save_pdb(f"{output_name}_{i}.pdb")

def batch_annotate_spacegroup(input_name, max_frame, sg):
    """
    Note: sg must be a string here!
    """
    for frame_id in range(max_frame):
        struct = gemmi.read_pdb(f"{input_name}_{frame_id}.pdb")
        struct.spacegroup_hm = sg
        struct.write_pdb(f"{input_name}_{frame_id}.pdb")

def batch_fmodel(input_name, max_frame, resolution=1.5,
                 phenix_command='source /usr/local/phenix-1.20.1-4487/phenix_env.sh; phenix.fmodel'):
    for frame_id in range(max_frame):
        ret = subprocess.run(phenix_command + f' {input_name}_{frame_id}.pdb high_resolution={resolution}')

def average_structure_factors(input_name):
    dataset = rs.read_mtz(f"{input_name}_0.pdb.mtz")
    n_reflections = dataset.shape[0]
    n_frames = 100

    complex_reflections = np.zeros(n_reflections, dtype='complex128')

    for frame in range(n_frames):
        dataset = rs.read_mtz(f"{input_name}_{frame}.pdb.mtz")
        complex_reflections = complex_reflections * (1 - 1/(frame + 1)) + np.array([amp*np.exp(np.pi*phase/180 * 1j) for [amp, phase] in dataset.to_numpy()]) / (frame + 1)


    dataset[:] = np.stack([np.abs(complex_reflections), np.angle(complex_reflections) / np.pi * 180]).T
    dataset.infer_mtz_dtypes(inplace = True)
    dataset.write_mtz(f"{input_name}_avg.mtz")

def generate_flexible_fun(f, reader, writer):
    """
    Transforms f(input)                                      -> retval
            to f_flex(generalized_input, generalized_output) -> None / retval'
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
                writer(f(reader(i)), output)
    return f_flex

def compute_net_dipole_moment(partial_charges, input=None, output=None):
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
    return (generate_flexible_fun(core, reader, writer))(input, output)

