#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e

yarn install

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
git stash

for tag in $(git tag); do
  git checkout $tag -- docs
  if [[ $(git ls-files | grep docs) ]]; then
    yarn run docusaurus docs:version $tag
  fi
done

git checkout $BRANCH
git stash pop

yarn build