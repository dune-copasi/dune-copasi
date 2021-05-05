#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e


# install cmake
curl -OL https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2-linux-x86_64.tar.gz
tar -xzf cmake-3.20.2-linux-x86_64.tar.gz
CMAKE=${PWD}/cmake-3.20.2-linux-x86_64/bin/cmake

yarn install

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
git stash

for tag in $(git tag); do
  echo "----------Checking out $tag-----------"
  if git checkout $tag -- docs; then
    echo "----------Inspecting files $tag-----------"
    if [[ $(git ls-files | grep docs) ]]; then
      git checkout $tag -- ../../dune
      mkdir docs/doxygen
      pushd docs/doxygen
      ${CMAKE} ../../../doxygen && make
      popd
      echo "----------Versioning $tag-----------"
      yarn run docusaurus docs:version $tag
      rm -rf docs/doxygen
    fi
  fi
done

echo "----------Get back to branch $BRANCH-----------"
git checkout $BRANCH -- docs
mkdir docs/doxygen
pushd docs/doxygen
${CMAKE} ../../../doxygen && make
popd

echo "----------Build yarn-----------"
yarn build
