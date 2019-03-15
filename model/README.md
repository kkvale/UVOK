# Prepare Metos3D UVOK model sources

Most of the UVOK model codes are part of the official [UVic](https://github.com/OSU-CEOAS-Schmittner/UVic2.9) model.
Others are copied from the official [TMM/UVOK repository](https://github.com/samarkhatiwala/tmm/tree/master/models/current/uvok1.0).

## Download UVic model sources

```sh
cd UVOK/uvic/
git clone https://github.com/OSU-CEOAS-Schmittner/UVic2.9.git
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz
tar xfz Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz
```

## Prepare UVic model files

```sh
echo '#include "../UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h
python prepare_uvic_files.py
```

## Clean up:

```sh
rm -fr include_uvok_tmm_options.h
rm -fr Kiel_Jan_2017_updates_to_UVIC2.9/
rm -fr Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz 
rm -fr UVic2.9/
```








https://github.com/samarkhatiwala/tmm.git
control.in
UVOK_TMM_OPTIONS.h

https://github.com/samarkhatiwala/tmm/blob/master/models/current/uvok1.0/runscripts/control.in
https://github.com/samarkhatiwala/tmm/blob/master/models/current/uvok1.0/src/UVOK_TMM_OPTIONS.h



## Prepare directory structure

```sh
mkdir UVOK
mkdir UVOK/uvic
```


## Download UVOK model sources

```sh
git clone https://github.com/samarkhatiwala/tmm.git
```

## Prepare UVOK model files

```sh
echo '#include "UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h

cp tmm/models/current/uvok1.0/runscripts/control.in UVOK/.
cp tmm/models/current/uvok1.0/src/UVOK_TMM_OPTIONS.h UVOK/.
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_ini.F > UVOK/uvok_ini.F
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_calc.F > UVOK/uvok_calc.F
```
