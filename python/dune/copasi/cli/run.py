import os
import sys
import subprocess

def run(args):
    if not os.path.isfile(args["config"]):
        raise IOError("Configuration file {} not found".format(args["config"]))
    try:
        subprocess.check_call([args["exec"], args["config"]])
    except subprocess.CalledProcessError:
        print("Error while running dune-copasi")
        sys.exit(1)