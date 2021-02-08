#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e

yarn install

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
git stash

for tag in $(git tag); do
  echo "----------Checking out $tag-----------"
  if git checkout $tag -- docs; then
    echo "----------Inspecting filder $tag-----------"
    if [[ $(git ls-files | grep docs) ]]; then
      echo "----------Versioning $tag-----------"
      yarn run docusaurus docs:version $tag
    fi
  fi
done

echo "----------Get back to branch $BRANCH-----------"
git checkout $BRANCH -- docs

echo "----------Build yarn-----------"
yarn build
