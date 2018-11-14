import mdtraj
import numpy as np

def rmsd(traj, ref, atom_indices, ref_atom_indices=None,
         align_atom_indices=None, align_ref_atom_indices=None,
         frameref=0):
    """
    RMSD observable with support for aligning on different sets of atoms 
    than used for RMSD calculation.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory for analysis. The RMSD is computed for all frames in
        this trajectory.
    ref : mdtraj.Trajectory
        Reference trajectory for computing RMSDs with respect to. This
        can be the same trajectory as traj, a different trajectory, or a
        PDB file that was loaded using mdtraj.
    atom_indices : np.ndarray
        Array of atomic indices in traj to use for RMSD calculation.
    ref_atom_indices : np.ndarray
        Array of atomic indices in ref to use for RMSD calculation. 
        Defaults to using indices provided as atom_indices. Must be the 
        same length as atom_indices
    align_atom_indices : np.ndarray
        Array of atomic indices in traj to use for alignment. If not
        provided, this defaults to atom_indices. 
    align_ref_atom_indices : np.ndarray
        Array of atomic indices in ref to use for alignment. If not
        provided, this defaults to ref_atom_indices. Must be the same 
        length as align_atom_indices
    frameref : int
        Frame in ref to compute RMSDs with respect to. Defaults to index
        0, which is either the first frame of a trajectory or the only
        frame for a PDB loaded by mdtraj.
    """
    # Check invariants
    if ref_atom_indices is not None:
        assert len(atom_indices) == len(ref_atom_indices)
    else:
        ref_atom_indices = atom_indices

    if align_ref_atom_indices is not None:
        assert align_atom_indices is not None
        assert len(align_atom_indices) == len(align_ref_atom_indices)
    else:
        if align_atom_indices is not None:
            align_ref_atom_indices = align_atom_indices
        else:
            align_atom_indices = atom_indices
            align_ref_atom_indices = ref_atom_indices

    # 1) Align to reference
    traj.superpose(ref, frameref, atom_indices=align_atom_indices,
                   ref_atom_indices=align_ref_atom_indices)

    # 2) Compute and return RMSD
    return np.sqrt(3*np.mean((traj.xyz[:, atom_indices, :] - ref.xyz[frameref, ref_atom_indices, :])**2, axis=(1,2)))
    
