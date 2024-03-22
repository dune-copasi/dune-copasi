#!/usr/bin/env python3

# This script updates version of the code in all places where it has to be hard coded
# and it updates the changelog. The version source of "truth" comes from dune.module in the root folder

import sys
import os
import re
import datetime
import argparse

# https://semver.org/spec/v2.0.0.html utility helpers
import semver
pattern = re.compile(r'^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$')

base_path = os.path.dirname(os.path.realpath(__file__))

def get_dune_module_version():
  dune_version_path = os.path.join(base_path,'../dune.module')
  with open(dune_version_path, "r") as version_file:
    for line in version_file:
      if line.startswith("Version:"):
        version = pattern.search(line.lstrip('Version: ').strip('\n'))
        assert version != None
        return semver.Version.parse(version.group())
  raise ValueError("No version in file " + dune_version_path)


def update_changelog(old_version, new_version):
  changelog_path = os.path.join(base_path,'../CHANGELOG.md')
  project_url = "https://gitlab.dune-project.org/copasi/dune-copasi"
  cl_new_version = new_version.finalize_version()
  cl_old_version = old_version.finalize_version()
  with open(changelog_path, "r") as changelog_file:
    new_changelog = ''
    for line in changelog_file.readlines():
      if "## [Unreleased]" in line:
        new_changelog += f"{line}\n"
        new_changelog += f"## [{cl_new_version}] ([git-diff][{cl_new_version}-diff]) - {datetime.date.today().isoformat()}\n"
      elif f"[Unreleased-diff]: {project_url}" in line:
        new_changelog += f"[Unreleased-diff]: {project_url}/compare/v{cl_new_version}...master\n"
        new_changelog += f"[{cl_new_version}-diff]: {project_url}/compare/v{cl_old_version}...v{cl_new_version}\n"
      elif f"[{cl_old_version}]: {project_url}" in line:
        new_changelog += f"[{cl_new_version}]: {project_url}/-/releases/v{cl_new_version}\n"
        new_changelog += line
      else:
        new_changelog += line

  with open(changelog_path, "w") as changelog_file:
    changelog_file.write(new_changelog)

  print(f"Changelog '{changelog_path}' was updated from '{cl_old_version}' to '{cl_new_version}', check manually that the update is correct")

def generic_version_update(path, old_version, new_version):
  with open(path, "r") as version_file:
    content = version_file.read()
  with open(path, "w") as version_file:
    version_file.write(re.sub(old_version, new_version, content))

def update_python_module(new_version, old_version):
  path = os.path.join(base_path,'../python/setup.py')
  prefix = 'version=\''
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  print(f"Python module '{path}' was updated from '{old_version}' to '{new_version}'")

def update_dune_module(old_version, new_version):
  path = os.path.join(base_path,'../dune.module')
  prefix = 'Version: '
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  print(f"Dune module '{path}' was updated from '{old_version}' to '{new_version}'")

def update_gitlab_ci(old_version, new_version):
  path = os.path.join(base_path,'../.gitlab-ci.yml')
  prefix = 'v'
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  print(f"GitLab CI config '{path}' was updated from '{old_version}' to '{new_version}'")

def update_npm(old_version, new_version):
  path = os.path.join(base_path,'../npm/package.json')
  prefix = '"version": "'
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  print(f"NPM package '{path}' was updated from '{old_version}' to '{new_version}'")

def update_docusaurus(old_version, new_version):
  path = os.path.join(base_path,'../doc/docusaurus/package.json')
  prefix = 'npm:@copasi/dune-copasi-wasm@'
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  prefix = '"version": "'
  generic_version_update(path, prefix+str(old_version), prefix+str(new_version))
  print(f"Docusaurus package '{path}' was updated from '{old_version}' to '{new_version}'")

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Version handler for dune-copasi", exit_on_error= False)
  parser.add_argument('-v', '--verbose', action='store_true')
  parser.add_argument("--override-version")
  parser.add_argument("--replace-major")
  parser.add_argument("--replace-minor")
  parser.add_argument("--replace-patch")
  parser.add_argument("--replace-prerelease")
  parser.add_argument("--replace-build")
  parser.add_argument("--finalize-version", nargs='?', const=True)
  parser.add_argument("--next", nargs='?', choices=['major', 'minor', 'patch', 'prerelease'])
  parser.add_argument("--update-dune-module", nargs='?', const=True)
  parser.add_argument("--update-changelog", nargs='?', const=True)
  parser.add_argument("--update-gitlab-ci", nargs='?', const=True)
  parser.add_argument("--update-npm", nargs='?', const=True)
  parser.add_argument("--update-docusaurus", nargs='?', const=True)
  # parser.add_argument("--update-python-module", nargs='?', const=True)
  parser.add_argument("--update-all", nargs='?', const=True)
  args = parser.parse_args()

  old_version = get_dune_module_version()
  if (args.verbose):
    print(f"Current version is '{old_version}'")

  new_version = old_version
  if args.override_version:
    new_version = semver.Version.parse(args.override_version)

  if args.next:
    new_version = new_version.next_version(args.next)

  if args.replace_major:
    new_version = new_version.replace(major=args.replace_major)
  if args.replace_minor:
    new_version = new_version.replace(minor=args.replace_minor)
  if args.replace_patch:
    new_version = new_version.replace(patch=args.replace_patch)
  if args.replace_prerelease:
    new_version = new_version.replace(prerelease=args.replace_prerelease)
  if args.replace_build:
    new_version = new_version.replace(build=args.replace_build)

  if args.finalize_version:
    new_version = new_version.finalize_version()

  if (args.verbose):
    print(f"New version is '{new_version}'\n")
  else:
    print(new_version)

  if args.update_changelog:
      update_changelog(old_version, new_version)
  if args.update_all or args.update_dune_module:
    update_dune_module(old_version, new_version)
  if args.update_all or args.update_gitlab_ci:
    update_gitlab_ci(old_version, new_version)
  if args.update_all or args.update_npm:
    update_npm(old_version, new_version)
  if args.update_all or args.update_docusaurus:
    update_docusaurus(old_version, new_version)
  # if args.update_all or args.update_python_module:
  #   update_python_module(old_version, new_version)

  if args.update_all and args.verbose:
    print("\nIf this is a definitive change, please update documentation and set a new tag in the repository\n")
    print(f"\tgit tag -s v{new_version}\n")
    print("\nAlso remember to forward the branch latest to the most recent tag\n")
    print(f"\tgit checkout latest && git merge --ff-only v{new_version}\n")
