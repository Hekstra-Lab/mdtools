# from setuptools import setup, find_packages
#
# # Get version number
# def getVersionNumber():
#     with open("mdtools/VERSION", 'r') as vfile:
#         version = vfile.read().strip()
#     return version
# __version__ = getVersionNumber()
#
# setup(
#     name='mdtools',
#     version=__version__,
#     author='Jack B. Greisman, Ziyuan Zhao',
#     author_email='greisman@g.harvard.edu, ziyuanzhao@fas.harvard.edu',
#     packages=find_packages(),
#     description='Tools for running MD simulation in OpenMM',
#     install_requires=[
#         "numpy >= 1.6"
#     ],
# )

import sys

sys.stderr.write(
    """
===============================
Unsupported installation method
===============================
matchmaps does not support installation with `python setup.py install`.
Please use `python -m pip install .` instead.
"""
)
sys.exit(1)


# The below code will never execute, however GitHub is particularly
# picky about where it finds Python packaging metadata.
# See: https://github.com/github/feedback/discussions/6456
#
# To be removed once GitHub catches up.

setup(  # noqa
    name="matchmaps",
    install_requires=[],
)