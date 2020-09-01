# Version number for mdtools
def getVersionNumber():
    import pkg_resources
    version = pkg_resources.require("mdtools")[0].version
    return version
__version__ = getVersionNumber()

def printVersions():
    """Print version numbers of mdtools and OpenMM"""
    from simtk import openmm
    print(f"mdtools: {__version__}")
    print(f"OpenMM: {openmm.__version__}")
    return

# Top-level API
from .prep.mdsystem import MDSystem
from .prep.solvatedmdsystem import SolvatedMDSystem
from .prep.latticemdsystem import  LatticeMDSystem
from .analysis.latticemdtrajectory import LatticeMDTrajectory
