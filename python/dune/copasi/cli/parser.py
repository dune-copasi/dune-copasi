import os
import argparse

try:
  from .cmake_variables import DUNE_COPASI_EXEC
except:
  pass

from .run import run
from .preprocess import preprocess_config

def get_parser():

  parser = argparse.ArgumentParser()

  subparsers = parser.add_subparsers(title="Commands", dest="command")
  subparsers.required = True

  parser_run = subparsers.add_parser('run',
      help="Run the dune-copasi main routine.",
      description="Run the dune-copasi main routine.")
  parser_run.add_argument('config',
                          help="Configuration file.")

  try:
    executable = os.environ['DUNE_COPASI_EXEC']
  except Exception as e:
    try:
      executable = DUNE_COPASI_EXEC
    except Exception as e:
      raise IOError("dune-copasi executable not found")

  parser_run.add_argument('--exec',
                           default=executable,
                           help="dune-copasi c++ executable.")
  
  parser_run.set_defaults(func=run)


  parser_preprocess = subparsers.add_parser('preprocess',
      help="Preprocess the config file.",
      description="Preprocess the config file.")
  parser_preprocess.add_argument('config',
                          help="Configuration file.")

  parser_preprocess.add_argument('new_config',
                          help="Preprocessed configuration file.")

  parser_preprocess.set_defaults(func=preprocess_config)

  return parser