# mdtools: Tools for Running MD Simulations in OpenMM

This package is intended to provide a useful framework for running MD
simulations using OpenMM. The goal of `mdtools` is to extend
OpenMM's funcitonality by facilitating common use cases with sensible default
parameter values. With `mdtools`, you can run an MD simulation of a protein in a
waterbox or the unit cell from a crystal structure with just a few lines of Python code.

For now, this repo is very much targeted to my current use cases, which is centered
around comparing MD simulations to crystallographic data.

# Installation

I have not yet made this package available on PyPI. Until then, `mdtools` can be
installed locally by cloning the repository and using the included `setup.py`:

```
git clone https://github.com/Hekstra-Lab/mdtools.git
cd mdtools
conda install -c omnia -c conda-forge openmm
conda install -c conda-forge mdtraj
python setup.py install
```

Note (from Ziyuan): I have added the requirements for `openmm` and `mdtraj` in 
`setup.py` but haven't tested it.


# Examples

Example scripts that make use of the provided libraries will be found in the
examples directory (Stay tuned).

# Changelog
* Added many utility functions in `mdtools.utils` for performing post-simulation
analyses on crystal MD systems.
* Support saving velocities and other data through additional reporters in 
`MDSystem.buildSimulation`.
* Support reporting force overflows at specific atoms in `MDSystem.equilibrate` 
and `MDSystem.simulate`, which are very helpful for debugging forcefield issues
during pilot MD simulation runs.
* Support bulk addition of specified amount of water during `squeeze` runs and 
set target volume to reach convergence more rapidly under prior knowledge about 
system density or volume.


