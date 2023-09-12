"""
trajectory.py: Provides LatticeMDTrajectory class extending MDTraj 
               Trajectory to provide custom analysis routines

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

import mdtraj
import numpy as np
import mdtools.utils as utils

# Problem, currently it doesn't seem possible to slice a LatticeMDTrajectory object if we do not override the slice method
class LatticeMDTrajectory(mdtraj.Trajectory):
    """
    LatticeMDTrajectory provides methods for the analysis of MD 
    simulations of crystal lattices.
    """
    def __init__(self, filename, sg=None):

        traj = mdtraj.load(filename)
        super().__init__(traj.xyz, traj.topology, traj.time,
                         traj.unitcell_lengths, traj.unitcell_angles)
        self.sg = sg

        # Validate unitcell information
        self._validateUnitcell()

    def _validateUnitcell(self):
        """
        Validate unitcell information is provided and do sanity checks
        """

        if not self._have_unitcell:
            raise AttributeError('Unitcell information is not provided')
        self._check_valid_unitcell()
        return

    def set_unit_cell_and_symmetry(unitcell_lengths=None, unitcell_angles=None, sg=None):
        if unitcell_lengths is not None:
            self.unitcell_lengths = unitcell_lengths
        if unitcell_angles is not None:
            self.unitcell_angles = unitcell_angles
        if sg is not None:
            self.sg = sg

    def revert_to_asu(self):
        """
        This method assumes that the trajectory is containing purely protein chains and then
        revert the coordinates of the atoms back to the ASU, according to the spacegroup
        represented by the trajectory.
        """
        if self.sg is None:
            raise Exception("You must specify a spacegroup when initializing the trajectory!")
        self.xyz = utils.trajectory_revert_to_asu(self, sg=self.sg, chain_ids=np.arange(self.n_chains)).xyz

    def smartWrapMolecule(self, indices):
        """
        This function applies periodic wrapping to a given set of atomic
        indices to prevent their center of mass from jumping by a unit
        cell length. Currently, it is assumed that the indices 
        correspond to a molecule -- meaning a set of atoms connected by
        bonds.

        Parameters
        ----------
        indices : list of ints
            Atomic indices of positions that should be wrapped together
        """

        # Compute geometric center of coordinates
        coms = self.xyz[:, indices, :].mean(axis=1)
    
        # Compute mask for integer unitcell adjustments
        mask = np.zeros(shape=(self.n_frames, 3))
    
        # X-axis
        x = self.unitcell_lengths[0, 0]
        mask[np.where(coms[:, 0] - coms[0, 0] < -1*x/2)[0], 0] = 1
        mask[np.where(coms[:, 0] - coms[0, 0] > x/2)[0], 0] = -1
    
        # Y-axis
        y = self.unitcell_lengths[0, 1]
        mask[np.where(coms[:, 1] - coms[0, 1] < -1*y/2)[0], 1] = 1
        mask[np.where(coms[:, 1] - coms[0, 1] > y/2)[0], 1] = -1
        
        # Z-axis
        z = self.unitcell_lengths[0, 2]
        mask[np.where(coms[:, 2] - coms[0, 2] < -1*z/2)[0], 2] = 1
        mask[np.where(coms[:, 2] - coms[0, 2] > z/2)[0], 2] = -1
    
        # Update trajectory coordinates
        self.xyz[:, indices, :] += (mask*self.unitcell_lengths).reshape(-1, 1, 3)    
    
        return 

    def smartWrapProtein(self):
        """
        Apply smart wrapping independently to each protein molecule in
        the MD system. For now, this method identifies proteins as 
        molecules with more than 100 atoms
        """

        for mol in self.topology.find_molecules():
            if len(mol) > 100:
                indices = [ atom.index for atom in mol ]
                self.smartWrapMolecule(indices)
                
