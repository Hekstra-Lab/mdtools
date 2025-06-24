# mdtools: Tools for Running MD Simulations in OpenMM

This package is intended to provide a useful framework for running MD
simulations using OpenMM. The goal of `mdtools` is to extend
OpenMM's funcitonality by facilitating common use cases with sensible default
parameter values. With `mdtools`, you can run an MD simulation of a protein in a
waterbox or the unit cell from a crystal structure with just a few lines of Python code.

For now, this repo is very much targeted to my current use cases, which is centered around comparing MD simulations to crystallographic data.

# Installation

## Pre-requisites to install via conda

The best course of action is to start off by installing `openmm`, `openmmforcefields`, and `openff-toolkit` via conda.

```shell
conda install -c omnia -c conda-forge openmm openmmforcefields openff-toolkit
```

### GPU specific installation

If you plan to use openmm/mdtools with GPU acceleration through CUDA (and you probably do!!) you'll need some extra goodies. As long as you run the above `conda install` command in the same GPU-containing environment that you plan to use for simulations, everything should be installed correctly automagically.

## pip-based installation of mdtools and other dependencies

I have not yet made this package available on PyPI. Until then, `mdtools` can be
installed via `git+` syntax:

```shell
pip install git+https://github.com/Hekstra-Lab/mdtools.git
```
or, you can clone the repo and install your local copy:

```
git clone https://github.com/Hekstra-Lab/mdtools.git
cd mdtools
pip install -e .
```


Note: building the package could be very slow in vanilla Conda due to the many dependencies that need to be resolved. Please consider installing [boa](https://github.com/mamba-org/boa) and use `conda mambabuild` for an immense speed-up. Also, use `mamba install` instead of `conda install` for additional speed-up.


# Examples

Example scripts that make use of the provided libraries will be found in the
examples directory (Stay tuned).

# Changelog
* Added and reformatted many docstrings for extensive documentation purpose.
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
* Several minor QoL changes.


