#!/usr/bin/env bash

# Loop over every tag and version documentation if available

set -e

# install cmake
curl -OL https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2-linux-x86_64.tar.gz
tar -xzf cmake-3.20.2-linux-x86_64.tar.gz
CMAKE=${PWD}/cmake-3.20.2-linux-x86_64/bin/cmake

yarn install

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')
echo "Branch: ${BRANCH}"

git stash

DOCUSAURUS_SOURCE=${PWD}
DOXYGEN_SOURCE=${PWD}/../doxygen
HEADERS_SOURCE=${PWD}/../../dune
DUNE_MODULE_FILE=${PWD}/../../dune.module

create_doxygen_html() {
  TAG=$1
  HTML_PATH=$2
  git checkout $TAG -- ${HEADERS_SOURCE} ${DUNE_MODULE_FILE}
  mkdir -p /tmp/doxygen && rm -rf /tmp/doxygen/*
  ${CMAKE} -S ${DOXYGEN_SOURCE} -B /tmp/doxygen
  make -C /tmp/doxygen
  mkdir -p ${HTML_PATH}
  mv /tmp/doxygen/html/* ${HTML_PATH}
  git checkout $BRANCH -- ${HEADERS_SOURCE} ${DUNE_MODULE_FILE}
}

for tag in $(git tag); do
  echo "----------Checking out $tag-----------"
  if git checkout $tag -- docs; then
    echo "----------Inspecting files $tag-----------"
    if [[ $(git ls-files | grep docs) ]]; then
      echo "----------Versioning $tag-----------"
      [[ "$BRANCH"=="master" ]] && create_doxygen_html $tag tmp_build/docs/$tag/doxygen
      yarn run docusaurus docs:version $tag
    fi
  fi
done

mkdir -p tmp_build/docs/doxygen
[[ "$BRANCH"=="master" ]] && mv tmp_build/docs/$tag/doxygen/* tmp_build/docs/doxygen

echo "----------Get back to branch $BRANCH-----------"
git checkout $BRANCH -- docs
create_doxygen_html $BRANCH tmp_build/docs/next/doxygen

echo "----------Build yarn-----------"
yarn build

rsync -av tmp_build/ build/
rm -rf tmp_build
