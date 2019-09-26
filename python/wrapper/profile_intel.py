#!/usr/bin/env python

if __name__ == "__main__":

  import sys
  import os
  import subprocess

  from dune.testtools.wrapper.argumentparser import get_args

  # Parse the given arguments
  args = get_args()
  ini = args["ini"]
  report_name = os.path.splitext(ini)[0] + "_r{:04d}"
  count = 0
  while (os.path.exists(report_name.format(count)) ):
    count+=1
  report_name = report_name.format(count)
  sys.exit(subprocess.call(['amplxe-cl', '-collect hotspots', '-r', report_name, '--', args["exec"], args["ini"]]))
