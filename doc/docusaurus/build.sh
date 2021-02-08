#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e

yarn install

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
git stash

for tag in $(git tag); do
  if git checkout $tag -- docs; then
    if [[ $(git ls-files | grep docs) ]]; then
      yarn run docusaurus docs:version $tag
    fi
  fi
done

git checkout $BRANCH
git stash pop &> /dev/null

yarn build