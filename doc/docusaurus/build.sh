#!/usr/bin/env bash

yarn install

set +e

BRANCH=$(cat ../../.git/HEAD | awk -F '/' '{print $NF}')

for tag in $(git tag)
  do yarn run docusaurus docs:version $tag
done

git checkout $BRANCH

yarn build