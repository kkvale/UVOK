# Metos3D UVOK

How to prepare the model file?

```sh
cd prepare_uvok_files/
git clone https://github.com/jpicau/UVic2.9.git
git clone https://github.com/samarkhatiwala/tmm.git
cd ../

python prepare_uvok_files.py

source:
uvok_calc.F............. prepare_uvok_files/tmm/models/current/uvok1.0/src/uvok_calc.F
gasbc.F................. prepare_uvok_files/UVic2.9/updates/03/source/common/gasbc.F
gosbc.F................. prepare_uvok_files/UVic2.9/updates/03/source/common/gosbc.F
setvbc.F................ prepare_uvok_files/UVic2.9/updates/03/source/mom/setvbc.F
tracer.F................ prepare_uvok_files/UVic2.9/updates/03/source/mom/tracer.F
uvok_stubs.F............ prepare_uvok_files/tmm/models/current/uvok1.0/src/uvok_stubs.F
uvok_diags_mod.F90...... prepare_uvok_files/tmm/models/current/uvok1.0/src/uvok_diags_mod.F90
co2calc.F............... prepare_uvok_files/UVic2.9/updates/02/source/common/co2calc.F
npzd_src.F.............. prepare_uvok_files/UVic2.9/updates/03/source/mom/npzd_src.F
uvok_ini.F.............. prepare_uvok_files/tmm/models/current/uvok1.0/src/uvok_ini.F
UVic_ESCM.F............. prepare_uvok_files/UVic2.9/updates/03/source/common/UVic_ESCM.F
iomngr.F................ prepare_uvok_files/UVic2.9/updates/03/source/common/iomngr.F
header:
UVOK_TMM_OPTIONS.h...... prepare_uvok_files/tmm/models/current/uvok1.0/src/UVOK_TMM_OPTIONS.h
size.h.................. prepare_uvok_files/UVic2.9/updates/03/source/common/size.h
param.h................. prepare_uvok_files/UVic2.9/source/common/param.h
pconst.h................ prepare_uvok_files/UVic2.9/source/common/pconst.h
stdunits.h.............. prepare_uvok_files/UVic2.9/source/common/stdunits.h
coord.h................. prepare_uvok_files/UVic2.9/source/common/coord.h
csbc.h.................. prepare_uvok_files/UVic2.9/updates/02/source/common/csbc.h
grdvar.h................ prepare_uvok_files/UVic2.9/source/common/grdvar.h
levind.h................ prepare_uvok_files/UVic2.9/source/common/levind.h
mw.h.................... prepare_uvok_files/UVic2.9/updates/03/source/mom/mw.h
scalar.h................ prepare_uvok_files/UVic2.9/source/common/scalar.h
tmngr.h................. prepare_uvok_files/UVic2.9/source/common/tmngr.h
npzd.h.................. prepare_uvok_files/UVic2.9/updates/03/source/common/npzd.h
calendar.h.............. prepare_uvok_files/UVic2.9/source/common/calendar.h
diaga.h................. prepare_uvok_files/UVic2.9/source/mom/diaga.h
ice.h................... prepare_uvok_files/UVic2.9/updates/02/source/common/ice.h
atm.h................... prepare_uvok_files/UVic2.9/updates/02/source/embm/atm.h
switch.h................ prepare_uvok_files/UVic2.9/updates/03/source/common/switch.h
cembm.h................. prepare_uvok_files/UVic2.9/updates/02/source/common/cembm.h
insolation.h............ prepare_uvok_files/UVic2.9/source/embm/insolation.h
solve.h................. prepare_uvok_files/UVic2.9/source/embm/solve.h
accel.h................. prepare_uvok_files/UVic2.9/source/common/accel.h
cregin.h................ prepare_uvok_files/UVic2.9/source/common/cregin.h
emode.h................. prepare_uvok_files/UVic2.9/source/common/emode.h
hmixc.h................. prepare_uvok_files/UVic2.9/updates/03/source/common/hmixc.h
timeavgs.h.............. prepare_uvok_files/UVic2.9/updates/03/source/mom/timeavgs.h
vmixc.h................. prepare_uvok_files/UVic2.9/source/common/vmixc.h
isopyc.h................ prepare_uvok_files/UVic2.9/updates/03/source/common/isopyc.h
fdift.h................. prepare_uvok_files/UVic2.9/source/mom/fdift.h
cnep.h.................. prepare_uvok_files/UVic2.9/source/common/cnep.h
cprnts.h................ prepare_uvok_files/UVic2.9/source/common/cprnts.h
diag.h.................. prepare_uvok_files/UVic2.9/updates/02/source/common/diag.h
fwa.h................... prepare_uvok_files/UVic2.9/source/common/fwa.h
iounit.h................ prepare_uvok_files/UVic2.9/source/common/iounit.h
mtlm.h.................. prepare_uvok_files/UVic2.9/updates/02/source/common/mtlm.h
sed.h................... prepare_uvok_files/UVic2.9/updates/02/source/common/sed.h
stab.h.................. prepare_uvok_files/UVic2.9/source/common/stab.h
veg.h................... prepare_uvok_files/UVic2.9/updates/02/source/embm/veg.h
iomngr.h................ prepare_uvok_files/UVic2.9/source/common/iomngr.h

```

