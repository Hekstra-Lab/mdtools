from setuptools import setup

setup(
    name='mdtools',
    version='0.1.0',
    author='Jack B. Greisman',
    author_email='greisman@g.harvard.edu',
    packages=['mdtools'],
    description='Tools for running MD simulation in OpenMM',
    install_requires=[
        "numpy >= 1.6",
        "mdtraj >= 1.9"
    ],
)
