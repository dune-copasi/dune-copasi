#!/bin/python3

# This script updates version of the code in all places where it has to be hard coded
# and it updates the changelog

# TODO update .gitlab.yml for dependency image
# TODO conditionally to update CHANGELOG

import sys
import os
import re
import datetime

# Regex to match and check version (taken from https://semver.org/spec/v2.0.0.html)
default_regex = r'^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$'

def check_version_format(version):
  if not re.match(default_regex,version):
    raise ValueError("\"" + version + "\" is not a valid version")

def get_version(path):
  with open(path, "r") as version_file:
    version = version_file.readline().rstrip("\n")
    version_file.close()
    check_version_format(version)
  return version

def update_version(path, old_version, new_version):
  with open(path, "r") as version_file:
    content = version_file.read()
  with open(path, "w") as version_file:
    version_file.write(re.sub(old_version, new_version, content))
  return 0

def main(argv):

  base_path = os.path.dirname(os.path.realpath(__file__))
  version_path = os.path.join(base_path,'../VERSION')
  python_version_path = os.path.join(base_path,'../python/setup.py')
  dune_version_path = os.path.join(base_path,'../dune.module')
  changelog_path = os.path.join(base_path,'../CHANGELOG.md')

  old_version = get_version(version_path)
  print("Old version is:", old_version)

  if len(argv) is 0:
    new_version = input("New version: ")
  elif len(argv) is 1:
    new_version = argv[0]
  else:
    raise IOError("Too many arguments!")

  check_version_format(new_version)

  # update basic version file
  update_version(version_path,old_version,new_version)

  # update python version file
  prefix = 'version=\''
  update_version(python_version_path,prefix+old_version,prefix+new_version)

  # update dune version file
  prefix = 'Version: '
  update_version(dune_version_path,prefix+old_version,prefix+new_version)

  project_url = "https://gitlab.dune-project.org/copasi/dune-copasi"

  with open(changelog_path, "r") as changelog_file:
    new_changelog = ''
    for line in changelog_file.readlines():
      if "## [Unreleased]" in line:
        new_changelog += line +"\n"
        new_changelog += "## [" + new_version+ "] ([git-diff][" + new_version + "-diff]) - " + datetime.date.today().isoformat() + "\n"
      elif "[Unreleased-diff]: " + project_url in line:
        new_changelog += "[Unreleased-diff]: " + project_url+"/compare/v"+new_version+"...master\n"
        new_changelog += "["+new_version+"-diff]: "+project_url+"/compare/v"+old_version+"...v"+new_version+"\n"
      elif "["+old_version+"]: "+project_url in line:
        new_changelog += "["+new_version+"]: "+project_url+"/-/releases/v"+new_version+"\n"
        new_changelog += line
      else:
        new_changelog += line

  with open(changelog_path, "w") as changelog_file:
    changelog_file.write(new_changelog)

  print("New version is:", new_version)

  print("Changelog was updated, check that the update is correct")

  print("If this is a definitive change, please update documentation and set a new tag in the repository\n")
  print("\tgit tag -s v"+new_version+"\n")
  print("\nAlso remember to forward the branch latest to the lastest tag\n")
  print("\tgit checkout latest && git merge --ff-only v"+new_version+"\n")

if __name__ == "__main__":
   main(sys.argv[1:])
