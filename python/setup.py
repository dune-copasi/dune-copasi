#!/usr/bin/env python

import sys
from setuptools import setup


def dune_copasi_scripts():
    return ['./scripts/dune-copasi']

setup(name='dune.copasi',
      version='0.0.0',
      namespace_packages=['dune'],
      description='Biochemical System Simulator COPASI in DUNE',
      author='Santiago Ospina De Los Ríos <santiago.ospina@iwr.uni-heidelberg.de>',
      author_email='no_mailinglist_yet@dune-copasi.de',
      url='http://copasi.org',
      packages=['dune.copasi',
                'dune.copasi.cli',
                ],
      install_requires=['sympy',
                        'configparser',
                        'argparse'],
      scripts=dune_copasi_scripts())
