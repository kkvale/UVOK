# Metos3D UVOK

## Prepare data

> **Note**
>
> We will use the [Miniconda Python ecosystem](https://docs.conda.io/en/latest/miniconda.html)
> to prepare the data.
>
> Install `conda` for Mac OS X, for instance:
>
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
# follow the instructions and let the installer modify your PATH variable
```
> Install `scipy` and `h5py` (HDF5):
>
```
conda install scipy
conda install h5py
```

### Download data and scripts:

```sh
cd data/
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVic_Kiel_increase_isopyc_diff.tar
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVic_Kiel_increase_isopyc_diff_model_data.tar
# curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/tmm_matlab_code.tar.gz
git clone https://github.com/samarkhatiwala/tmm.git
```

### Extract archives:

```sh
tar xf UVic_Kiel_increase_isopyc_diff.tar 
tar xf UVic_Kiel_increase_isopyc_diff_model_data.tar 
#mkdir tmm_matlab_code
#tar xfz tmm_matlab_code.tar.gz -C tmm_matlab_code
```

### Prepare directory structure:

```sh
mkdir geometry/
mkdir ini/
mkdir transport/
mkdir -p forcing/domain/
mkdir -p forcing/boundary/
```

### Create script `prepare_uvok_data.py` ...

### Run script:

```sh
python prepare_uvok_data.py 
```

### Clean up:

```sh
rm -fr UVic_Kiel_increase_isopyc_diff.tar UVic_Kiel_increase_isopyc_diff/
rm -fr UVic_Kiel_increase_isopyc_diff_model_data.tar UVic_Kiel_increase_isopyc_diff_model_data/
#rm -fr tmm_matlab_code.tar.gz tmm_matlab_code
rm -fr tmm
```

## Prepare model sources

### Download model sources

```sh
cd model/
git clone https://github.com/samarkhatiwala/tmm.git
```

### Prepare directory structure

```sh
cd model/
mkdir UVOK
mkdir UVOK/option
mkdir UVOK/uvic
```
### Prepare UVOK model files

```sh
echo '#include "UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h

cp 
cp tmm/models/current/uvok1.0/src/UVOK_TMM_OPTIONS.h UVOK/.
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_ini.F > UVOK/uvok_ini.F
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_calc.F > UVOK/uvok_calc.F
#cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_stubs.F > UVOK/uvok_stubs.F
#cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_diags_mod.F90 > UVOK/uvok_diags_mod.F90  
```

### Prepare UVic model files

```sh
cd UVOK/uvic/
echo '#include "../UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h

git clone https://github.com/OSU-CEOAS-Schmittner/UVic2.9.git
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz

tar xfz Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz

python prepare_uvic_files.py
```

### Clean up:

```sh
rm -fr tmm
rm -fr include_uvok_tmm_options.h

# we leave source directories in `uvic/` intact, to which the header files are linked  
#   UVic2.9/
#   Kiel_Jan_2017_updates_to_UVIC2.9/
```


