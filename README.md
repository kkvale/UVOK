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
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/tmm_matlab_code.tar.gz
git clone https://github.com/samarkhatiwala/tmm.git tmm_github
```

### Extract archives:

```sh
tar xf UVic_Kiel_increase_isopyc_diff.tar 
tar xf UVic_Kiel_increase_isopyc_diff_model_data.tar 
mkdir tmm_matlab_code
tar xfz tmm_matlab_code.tar.gz -C tmm_matlab_code
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

```
rm -fr UVic_Kiel_increase_isopyc_diff.tar UVic_Kiel_increase_isopyc_diff/
rm -fr UVic_Kiel_increase_isopyc_diff_model_data.tar UVic_Kiel_increase_isopyc_diff_model_data/
rm -fr tmm_matlab_code.tar.gz tmm_matlab_code
rm -fr tmm
```

## Prepare model sources


