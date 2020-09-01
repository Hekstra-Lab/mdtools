# Version number for mdtools
def getVersionNumber():
    import pkg_resources
    version = pkg_resources.require("mdtools")[0].version
    return version
__version__ = getVersionNumber()

from .prep.mdsystem import MDSystem
from .prep.solvatedmdsystem import SolvatedMDSystem
from .prep.latticemdsystem import  LatticeMDSystem
from .analysis.latticemdtrajectory import LatticeMDTrajectory
