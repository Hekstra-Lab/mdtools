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
git clone https://github.com/JBGreisman/mdtools.git
cd mdtools
conda install -c omnia -c conda-forge openmm
python setup.py install
```

# Examples

Example scripts that make use of the provided libraries will be found in the
examples directory (Stay tuned).


