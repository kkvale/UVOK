# Metos3D UVOK

This repository consists of the required resources to run the
[UVOK](https://github.com/samarkhatiwala/tmm/tree/master/models/current/uvok1.0) model with
[Metos3D](https://github.com/metos3d). See the particular subdirectories for details:

- [`data/`](data/) How to prepare the data?
- [`model/`](model/) How to prepare the model sources?
- [`test/`](test/) How to run Metos3D with UVOK?

> **Note:**
> We use the [Miniconda Python ecosystem](https://docs.conda.io/en/latest/miniconda.html)
> to prepare the data. In particular, we need `scipy` to access *older* and
> `h5py` to access *newer* Matlab files. Moreover, we use `numpy` arrays internally.
> However, the `numpy` module is installed along with `scipy`.

Install `conda` including **Python 3** for **Mac OS X**:

```sh
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
# follow the instructions and let the installer modify your PATH variable
```

Install `conda` including **Python 3** for **Linux**:

```sh
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# follow the instructions and let the installer modify your PATH variable
```

If `conda` has been installed, get `scipy` and `h5py` (HDF5):

```sh
conda install scipy
conda install h5py
```


