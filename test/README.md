# Metos3D UVOK test run

```
mistral.dkrz.de
rzcluster.rz.uni-kiel.de
nesh-fe.rz.uni-kiel.de
glogin.hlrn.de
```


PETSc is installed, tested PETSc 3.7 - 3.11

```sh
mkdir petsc
cd petsc
curl -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.11.0.tar.gz
tar xzf petsc-lite-3.11.0.tar.gz
cd petsc-3.11.0/
```

DKRZ:

```sh
module load intel/18.0.4
module load intelmpi/2018.4.274
./configure --with-debugging=0 --CC='mpicc -cc=icc' --CXX='mpicxx -cxx=icpc' --FC='mpiifort -fc=ifort' --with-blaslapack-dir=/sw/rhel6-x64/intel/intel-18.0.4/mkl
make PETSC_DIR=/mnt/lustre01/work/bk0993/b380415/playground/uvok/uvok/test/petsc/petsc-3.11.0 PETSC_ARCH=arch-linux-c-opt all
```

`petsc.env.sh`:

```
module load intel/18.0.4
module load intelmpi/2018.4.274
export PETSC_DIR=/mnt/lustre01/work/bk0993/b380415/playground/uvok/uvok/test/petsc/petsc-3.11.0
export PETSC_ARCH=arch-linux-c-opt
```

Metos3D is installed

```
python -c "$(curl -fsSL https://github.com/metos3d/metos3d/raw/master/install.py)"
```

> **Note:** The UVic model code is old. Back then no one seemed to care about 
> ...
> ... not to mention good scientific programming.
> ... clear interfaces, librarization
> ... `real` has no explicit type such as `real(4)` or `real(8)`,
> this has to be configured via a compiler switch,
> unfortunately, every compiler has its own switch for that, for instance,
> GNU, `gfortran`, `-fdefault-real-8`
> Intel, `ifort`, `-r8`
> IBM, `xlf`, `-qrealsize=8`
> and so on and so forth ...
> even worse, if using MPI compiler wrappers (along wit a specific MPI library),
> for instance, Open MPI, MPICH, MVAPICH, ...   
> they are not distiguishable from the shell command name, and 
> thus, you have to swallow the pill, you have to provide them by yourself, 

```sh
ln -s ~/.metos3d/metos3d/Makefile 
ln -s ~/.metos3d/simpack/
ln -s ../model/
ln -s ../model/UVOK/control.in 
ln -s ../data/
mkdir work/

. petsc/petsc.env.sh
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

openmpi, does not work, 

```
python2.7 ./configure --with-x=0 --CC=mpicc --CXX=mpic++ --FC=mpif90 --FFLAGS='-Wl,-rpath,/Users/jpicau/miniconda3/lib' --with-debugging=0
 make PETSC_DIR=/Users/jpicau/Documents/development/uvok/test/petsc/petsc-3.10.4 PETSC_ARCH=arch-darwin-c-opt all
```

dowload mpich, only one that works until now 

```
python2.7 ./configure --with-x=0 --with-debugging=0 --CC=gcc --CXX=g++ --FC=gfortran --FFLAGS='-Wl,-rpath,/Users/jpicau/miniconda3/lib' --download-mpich=1
make PETSC_DIR=/Users/jpicau/Documents/development/uvok/test/petsc/petsc-3.10.4 PETSC_ARCH=arch-darwin-c-opt all
```
