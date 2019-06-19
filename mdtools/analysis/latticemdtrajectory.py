"""
trajectory.py: Provides LatticeMDTrajectory class extending MDTraj 
               Trajectory to provide custom analysis routines

Author: Jack Greisman <greisman@g.harvard.edu>
"""
__author__  = "Jack Greisman"
__version__ = "1.0"

import mdtraj

class LatticeMDTrajectory(mdtraj.Trajectory):
    """
    LatticeMDTrajectory provides methods for the analysis of MD 
    simulations of crystal lattices.
    """
    def __init__(self, filename):
        traj = mdtraj.load(filename)
        super().__init__(traj.xyz, traj.topology, traj.time,
                         traj.unitcell_lengths, traj.unitcell_angles)

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
    
