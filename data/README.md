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
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVicKielIncrIsopycDiff.tar
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/UVicKielIncrIsopycDiff_model_data.tar
git clone https://github.com/samarkhatiwala/tmm.git
```

## Extract archives:

```sh
tar xf UVicKielIncrIsopycDiff.tar 
tar xf UVicKielIncrIsopycDiff_model_data.tar 
```

## Run script:

```sh
python prepare_uvok_data.py 
```

## Clean up:

```sh
rm -fr UVicKielIncrIsopycDiff.tar UVicKielIncrIsopycDiff/
rm -fr UVicKielIncrIsopycDiff_model_data.tar UVicKielIncrIsopycDiff_model_data/
rm -fr tmm
```


