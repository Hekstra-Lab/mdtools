# mdtools: Tools for Running MD Simulations in OpenMM

This package is intended to provide a useful framework for running MD
simulations using OpenMM. The goal of `mdtools` is to extend
OpenMM's funcitonality by facilitating common use cases with sensible default
parameter values. With `mdtools`, you can run an MD simulation of a protein in a
waterbox or the unit cell from a crystal structure with just a few lines of Python code.

For now, this repo is very much targeted to my current use cases, which is centered around comparing MD simulations to crystallographic data.

# Installation

I have not yet made this package available on PyPI. Until then, `mdtools` can be
installed locally by cloning the repository and using the included `setup.py`:

```
git clone https://github.com/Hekstra-Lab/mdtools.git
cd mdtools
python setup.py install
```

However, one would need to pre-install the following packages as specified in `conda-recipe/meta.yaml`:

```
- openmm>=7
- reciprocalspaceship=1.0.0 
- gemmi 
- mdtraj
- matplotlib
```

For example, openmm and mdtraj could be installed from conda:
```
conda install -c omnia -c conda-forge openmm
conda install -c conda-forge mdtraj
```

Alternatively, build the package by yourself and install all the above dependencies automatically with Conda:
```
conda create -n NEW_CONDA_ENV_NAME python=DESIRED_PYTHON_VERSION # >=3.7 
conda activate NEW_CONDA_ENV_NAME
conda build . --python=DESIRED_PYTHON_VERSION #choose from 3.7, 3.8, 3.9, 3.10, 3.11
conda install --use-local mdtools
```

Note: building the package could be very slow in vanilla Conda due to the many dependencies that need to be resolved. Please consider installing [boa](https://github.com/mamba-org/boa) and use `conda mambabuild` for an immense speed-up. Also, use `mamba install` instead of `conda install` for additional speed-up.

The above Conda build and install method has been tested on arm64 OSX and also on FASRC clusters. 


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


