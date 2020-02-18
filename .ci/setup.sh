# dependencies setup script for Travis and AppVeyor CI

DUNE_VERSION="2.7"

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNECONTROL: ${DUNECONTROL}"
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
  git clone -b support/dune-copasi --depth 1 --recursive https://gitlab.dune-project.org/copasi/$repo.git
done

# python virtual environment does not work in windows yet
if [[ ! $MSYSTEM ]]; then
	git clone -b feature/allow-multidomain-vtk-compare-to-have-same-thresholds https://gitlab.dune-project.org/quality/dune-testtools.git
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
wget https://gist.githubusercontent.com/lkeegan/059984b71f8aeb0bbc062e85ad7ee377/raw/e9c7af42c47fe765547e60833a72b5ff1e78123c/cmake-patch.txt
echo '' >> cmake-patch.txt
git apply cmake-patch.txt
# another patch for missing header in cmake install list
git apply ../dune-copasi/.ci/dune-common.patch
cd ../

cd dune-logging
git apply ../dune-copasi/.ci/dune-logging.patch
cd ../

ls

# python virtual environment does not work in windows yet
if [[ ! $MSYSTEM ]]; then
	${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --module=dune-testtools all
fi

for repo in dune-testtools dune-logging dune-pdelab dune-multidomaingrid
do
  ${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --module=$repo all
done