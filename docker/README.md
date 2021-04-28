## Dockerfiles for dune-copasi

There is one docker file which is of most importance:
(`dune-copasi.dockerfile`)[dune-copasi.dockerfile]. This file is able to
perform all the stages required to produce the `dune-copasi` executables
based on the local source files. It includes all the stages described in
[CI](../.ci) but also installs and test the generated packages into a fresh
encironment. Its purpose is two-fold:

* Continuous Integration: Its intermediate target `setup-env` is the starting point of the GitLab CI jobs.
* Development: Build changes producedd on local environments without need of installing `dune-copasi` dependencies.

## Development with dockerfiles

You only need `git` and `docker`, it cannot be simpler!

```
# clone main repository
git clone ssh://git@gitlab.dune-project.org:22022/copasi/dune-copasi.git

# explore source files and make local changes
cd dune-copasi

# build a new container with the local changes
docker build -t dune-copasi .

# test changes with your local config files (be sure that $PWD gives docker read/write acces)
docker run -v $PWD dune-copasi config.ini
```

## Production containers

The file [`deploy.dockerfile`](deploy.dockerfile) basically defines the last stage
of the development container. The CI passes its packaged results and this container
installs and tests the the packages. Finally, it is deployied to our registry. Its
main purpose is to create docker containers with executables usable by end-users.

## Jupiter Notebooks

The file [`jupyter.dockerfile`](jupyter.dockerfile) sets up an environment where
jupyter notebooks with C++ kernels and `dune-copasi` are available. This is
useful in order to explore the C++ API without the need of installing the
complete stack.

```
# clone main repository
git clone ssh://git@gitlab.dune-project.org:22022/copasi/dune-copasi.git

# set current path to dune-copasi root directory
cd dune-copasi

# build a local container for jupyter notebooks
docker build -t dune-copasi-jupiter -f docker/jupyter.dockerfile .

# run container, publish its contents to 8888 port, and follow instructions in the terminal to open the notebooks
docker run -ti -p 8888:8888 dune-copasi-jupiter
```
