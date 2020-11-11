#!/usr/bin/env python

import sys
from setuptools import setup

def dune_copasi_scripts():
    return ['./wrapper/profile_intel.py']

setup(name='dune.copasi',
  version='0.4.0-git',
  namespace_packages=['dune'],
  description='DUNE Biochemical System Simulator for COPASI',
  author='Santiago Ospina De Los Ríos',
  author_email='santiago.ospina@iwr.uni-heidelberg.de',
  url='http://copasi.org',
  packages=['dune.copasi'],
  install_requires=['dune.testtools'],
  scripts=dune_copasi_scripts()
)