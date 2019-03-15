# Prepare Metos3D UVOK data  

## Prepare directory structure:

```sh
mkdir geometry/
mkdir ini/
mkdir transport/
mkdir -p forcing/domain/
mkdir -p forcing/boundary/
```

## Download data and scripts:

```sh
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVic_Kiel_increase_isopyc_diff.tar
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVic_Kiel_increase_isopyc_diff_model_data.tar
git clone https://github.com/samarkhatiwala/tmm.git
```

## Extract archives:

```sh
tar xf UVic_Kiel_increase_isopyc_diff.tar 
tar xf UVic_Kiel_increase_isopyc_diff_model_data.tar 
```

## Run script:

```sh
python prepare_uvok_data.py 
```

## Clean up:

```sh
rm -fr UVic_Kiel_increase_isopyc_diff.tar UVic_Kiel_increase_isopyc_diff/
rm -fr UVic_Kiel_increase_isopyc_diff_model_data.tar UVic_Kiel_increase_isopyc_diff_model_data/
rm -fr tmm
```


