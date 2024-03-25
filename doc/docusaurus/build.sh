#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e

# install cmake
curl -OL https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2-linux-x86_64.tar.gz
tar -xzf cmake-3.20.2-linux-x86_64.tar.gz
CMAKE=${PWD}/cmake-3.20.2-linux-x86_64/bin/cmake

git lfs install
git lfs pull

# install dependencies from yarn lock
yarn install

# upgrade dune-copasi-wasm to latest uploaded version
yarn upgrade dune-copasi-wasm-git

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
echo "Branch: ${BRANCH}"

git stash

echo "----------Build yarn-----------"
yarn build
