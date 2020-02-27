from setuptools import setup, find_packages

setup(
    name='mdtools',
    version='0.1.1',
    author='Jack B. Greisman',
    author_email='greisman@g.harvard.edu',
    packages=find_packages(),
    description='Tools for running MD simulation in OpenMM',
    install_requires=[
        "numpy >= 1.6",
        "mdtraj >= 1.9",
        "openmm",
    ],
)
