# Metos3D UVOK test run


PETSc is installed
Metos3D is installed


```sh
ln -s ~/.metos3d/metos3d/Makefile 
ln -s ~/.metos3d/simpack/
ln -s ../model/
ln -s ../data/
mkdir work/

make BGC=model/UVOK clean
make BGC=model/UVOK

# ifort, FFLAGS+=-r8
# gfortran, FFLAGS+=-fdefault-real-8
# xlf, FPPFLAG+=-WF, FFLAGS+=-qrealsize=8 -qzerosize


```

PETSc

```sh
mkdir petsc
cd petsc/
curl -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.10.4.tar.gz
touch petsc.env.sh

tar xfz petsc-lite-3.10.4.tar.gz

cd petsc-3.10.4/

# debug
python2.7 ./configure --with-mpi=0 --with-x=0 --FC=gfortran --FFLAGS='-Wl,-rpath,/Users/jpicau/miniconda3/lib'
make PETSC_DIR=/Users/jpicau/Documents/development/uvok/test/petsc/petsc-3.10.4 PETSC_ARCH=arch-darwin-c-debug all

# opt
python2.7 ./configure --with-mpi=0 --with-x=0 --FC=gfortran --FFLAGS='-Wl,-rpath,/Users/jpicau/miniconda3/lib' --with-debugging=0
make PETSC_DIR=/Users/jpicau/Documents/development/uvok/test/petsc/petsc-3.10.4 PETSC_ARCH=arch-darwin-c-opt all

```


