from setuptools import setup, find_packages

# Get version number
def getVersionNumber():
    with open("mdtools/VERSION", 'r') as vfile:
        version = vfile.read().strip()
    return version
__version__ = getVersionNumber()

setup(
    name='mdtools',
    version=__version__,
    author='Jack B. Greisman',
    author_email='greisman@g.harvard.edu',
    packages=find_packages(),
    description='Tools for running MD simulation in OpenMM',
    install_requires=[
        "numpy >= 1.6",
    ],
)
