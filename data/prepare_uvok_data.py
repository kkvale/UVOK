#
# Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
# Copyright (C) 2019  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import scipy.sparse as spsp
import scipy.io as spio
import h5py as h5
import numpy as np
import numpy.ma as ma

### -----------------------------------------------------------------------------------------------

def writePETScVector(filename, vec):
    with open(filename, 'wb') as f:
        # vec id, nvec, data
        nvec, = vec.shape
        np.array([1211214, nvec], dtype='>i4').tofile(f)
        np.array(vec, dtype='>f8').tofile(f)

def writePETScMatrix(filename, mat):
    with open(filename, 'wb') as f:
        # mat id, nrow, ncol, nnz
        nrow, ncol = mat.shape
        nnz = mat.nnz
        np.array([1211216, nrow, ncol, nnz], dtype='>i4').tofile(f)
        # nnzrow, colind, data
        np.array(np.diff(mat.indptr), dtype='>i4').tofile(f)
        np.array(mat.indices, dtype='>i4').tofile(f)
        np.array(mat.data, dtype='>f8').tofile(f)

### -----------------------------------------------------------------------------------------------

# grid.mat
# file ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat
# ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat: Hierarchical Data Format (version 5) with 512 bytes user block
with h5.File('./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat') as f:
    
    # mask
    msk = (f['bathy'][...] != 1.)
    nz, ny, nx = msk.shape

#    # geometry/landSeaMask.petsc
#    lsm = spsp.csr_matrix(f['ideep'][...])
#    writePETScMatrix('geometry/landSeaMask.petsc', lsm)
#
#    # geometry/volumes.petsc
#    vol = ma.array(f['dv'][...], mask=msk)
#    vol = np.reshape(vol, (nz, ny*nx)).transpose()
#    vol = np.reshape(vol, (ny, nx, nz))
#    vol = vol.compressed()
#    writePETScVector('geometry/volumes.petsc', vol)
#
#    # forcing/domain/dz.petsc
#    dz = ma.array(f['dz'][...], mask=msk)
#    dz = np.reshape(dz, (nz, ny*nx)).transpose()
#    dz = np.reshape(dz, (ny, nx, nz))
#    dz = dz.compressed()
#    # scale, m to cm
#    dz = 100*dz
#    writePETScVector('forcing/domain/dz.petsc', dz)

    # forcing/boundary/latitude.petsc
    y = f['y'][...].transpose()
    phi = np.zeros((ny, nx))
    phi[...] = y[:]
    phi = ma.array(phi, mask=msk[0,:,:])
    phi = phi.compressed()
    writePETScVector('forcing/boundary/latitude.petsc', phi)

## init/*ini.petsc
#names = ['dic','c14','alk','o2','po4','phyt','zoop','detr','no3','diaz']
#for name in names:
#    filenamein = 'tmm_github/models/current/uvok1.0/matlab/InitialConditionProfiles/' + name + '.dat'
#    # vec1d
#    vec1d = np.loadtxt(filenamein)[:,1]
#    # vec3d
#    vec3d = np.zeros(msk.shape)
#    vec3d[...] = vec1d[:, np.newaxis, np.newaxis]
#    # mask, reshape
#    vec = ma.array(vec3d, mask=msk)
#    vec = np.reshape(vec, (nz, ny*nx)).transpose()
#    vec = np.reshape(vec, (ny, nx, nz))
#    vec = vec.compressed()
#    # write
#    filenameout = 'ini/' + name + 'ini.petsc'
#    writePETScVector(filenameout, vec)








#./UVic_Kiel_increase_isopyc_diff_model_data/GCM/basin_mask.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/GCM/Salt_gcm.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/GCM/FreshWaterForcing_gcm.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/GCM/Theta_gcm.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/config_data.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/BiogeochemData/ice_fraction.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/BiogeochemData/wind_speed.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff_model_data/BiogeochemData/UVOK_input_data.mat: Matlab v5 mat-file (little endian) version 0x0100
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_annualmean.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_08.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_09.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_07.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_12.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_06.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_10.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_04.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_05.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_11.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_01.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_02.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_03.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/boxes.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/tracer_tiles.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/basis_functions.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/matrix_extraction_run_data.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/boxnum.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/links.mat: Hierarchical Data Format (version 5) with 512 bytes user block
#./UVic_Kiel_increase_isopyc_diff/Matrix1/Data/profile_data.mat: Hierarchical Data Format (version 5) with 512 bytes user block

