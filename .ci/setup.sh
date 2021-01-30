#!/usr/bin/env bash

# dependencies setup script for Travis and AppVeyor CI
# requisites
#   * DUNE_OPTIONS_FILE is defined

set -e

DUNE_VERSION="2.7"

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNE_OPTIONS_FILE: ${DUNE_OPTIONS_FILE}"
cat ${DUNE_OPTIONS_FILE}
echo "PWD: $PWD"

which g++
g++ --version

which gcc
gcc --version

which python
python --version

which python3
python3 --version

which pip
pip --version

which pip3
pip3 --version

which cmake
cmake --version

# download Dune dependencies
for repo in core/dune-common core/dune-geometry core/dune-grid core/dune-istl core/dune-localfunctions staging/dune-functions staging/dune-uggrid
do
  git clone -b releases/$DUNE_VERSION --depth 1 --recursive https://gitlab.dune-project.org/$repo.git
done
for repo in dune-logging dune-typetree dune-pdelab dune-multidomaingrid
do
  git clone -b support/dune-copasi-v0.3 --depth 1 --recursive https://gitlab.dune-project.org/copasi/$repo.git
done

# python virtual environment does not work in windows yet
if [[ ! $MSYSTEM ]]; then
	git clone -b releases/$DUNE_VERSION https://gitlab.dune-project.org/quality/dune-testtools.git
fi

# on windows, symlinks from git repos don't work
# msys git replaces symlinks with a text file containing the linked file location
# so here we identify all such files, and replace them with the linked file
# note msys defines MSYSTEM variable: use this to check if we are on msys/windows
if [[ $MSYSTEM ]]; then
	rootdir=$(pwd)
		for repo in $(ls -d dune-*/)
		do
			echo "repo: $repo"
			cd $rootdir/$repo
			for f in $(git ls-files -s | awk '/120000/{print $4}')
			do
				dname=$(dirname "$f")
				fname=$(basename "$f")
				relf=$(cat $f)
				src="$rootdir/$repo/$dname/$relf"
				dst="$rootdir/$repo/$dname/$fname"
				echo "  - copying $src --> $dst"
				cp $src $dst
			done
		done
	cd $rootdir
fi

# patch cmake macro to avoid build failure when fortran compiler not found, e.g. on osx
cd dune-common
if [[ "$OSTYPE" == "darwin"* ]]; then
	wget https://gist.githubusercontent.com/lkeegan/059984b71f8aeb0bbc062e85ad7ee377/raw/e9c7af42c47fe765547e60833a72b5ff1e78123c/cmake-patch.txt
	echo '' >> cmake-patch.txt
	git apply cmake-patch.txt
fi

# another patch for missing header in cmake install list
git apply ../dune-copasi/.ci/dune-common.patch
cd ../


cd dune-logging
git apply ../dune-copasi/.ci/dune-logging.patch
cd ../

# hardcoded *ordered* dependencies!
MODULES="logging uggrid geometry grid localfunctions istl typetree functions pdelab multidomaingrid "
if [[ ! $MSYSTEM ]]; then
	MODULES+="testtools"
fi

./dune-common/bin/dunecontrol --opts=${DUNE_OPTIONS_FILE} --only=dune-common all
cmake --build dune-common/build-cmake/ --target install
rm -rf dune-common

# get cmake flags
CMAKE_FLAGS=
if test "x$DUNE_OPTIONS_FILE" != "x"; then
  CMAKE_FLAGS="$(. $DUNE_OPTIONS_FILE; eval echo \$CMAKE_FLAGS)"
fi

# get install prefix
for flag in $CMAKE_FLAGS; do
  [[ ${flag#-D} == CMAKE_INSTALL_PREFIX* ]] && CMAKE_INSTALL_PREFIX="${flag#-DCMAKE_INSTALL_PREFIX:PATH=}"
done

# define DUNECONTROL path
if test "x$DUNECONTROL" == "x"; then
	if test "x$CMAKE_INSTALL_PREFIX" != "x"; then
  	DUNECONTROL="$CMAKE_INSTALL_PREFIX/bin/dunecontrol"
	else
		DUNECONTROL=dunecontrol
	fi
fi

if ! command -v $DUNECONTROL &> /dev/null; then
	echo "ERROR: '$DUNECONTROL' cannot be executed"
	exit 1
fi

for module in $MODULES
do
	${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --only=dune-$module all
	cmake --build dune-$module/build-cmake/ --target install
	rm -rf dune-$module
done
