# Prepare Metos3D UVOK model sources

## Download UVOK model sources

```sh
git clone https://github.com/samarkhatiwala/tmm.git
```

## Prepare directory structure

```sh
mkdir UVOK
mkdir UVOK/uvic
```

## Prepare UVOK model files

```sh
echo '#include "UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h

cp tmm/models/current/uvok1.0/runscripts/control.in UVOK/.
cp tmm/models/current/uvok1.0/src/UVOK_TMM_OPTIONS.h UVOK/.
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_ini.F > UVOK/uvok_ini.F
cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_calc.F > UVOK/uvok_calc.F
#cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_stubs.F > UVOK/uvok_stubs.F
#cat include_uvok_tmm_options.h tmm/models/current/uvok1.0/src/uvok_diags_mod.F90 > UVOK/uvok_diags_mod.F90  
```

## Prepare UVic model files

```sh
cd UVOK/uvic/
echo '#include "../UVOK_TMM_OPTIONS.h"' > include_uvok_tmm_options.h

git clone https://github.com/OSU-CEOAS-Schmittner/UVic2.9.git
curl -O http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz

tar xfz Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz

python prepare_uvic_files.py
```

## Clean up:

```sh
rm -fr tmm
rm -fr include_uvok_tmm_options.h

# we leave source directories in `uvic/` intact, to which the header files are linked  
#   UVic2.9/
#   Kiel_Jan_2017_updates_to_UVIC2.9/
```



