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

def writePETScVector(filename, vec):
    with open(filename, 'wb') as f:
        # vec id, nvec
        nvec, = vec.shape
        np.array([1211214, nvec], dtype='>i4').tofile(f)
        # data
        np.array(vec, dtype='>f8').tofile(f)

def writePETScMatrix(filename, mat):
    with open(filename, 'wb') as f:
        # mat id, nrow, ncol, nnz
        nrow, ncol = mat.shape
        nnz = mat.nnz
        np.array([1211216, nrow, ncol, nnz], dtype='>i4').tofile(f)
        # nnzrow
        # colind
        # data
        np.array(np.diff(mat.indptr), dtype='>i4').tofile(f)
        np.array(mat.indices, dtype='>i4').tofile(f)
        np.array(mat.data, dtype='>f8').tofile(f)

# ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat, HDF5
with h5.File('./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat') as f:
    
    # landSeaMask.petsc
    lsm = spsp.csr_matrix(f['ideep'][...])
    writePETScMatrix('geometry/landSeaMask.petsc', lsm)
    
    # mask
    msk = (f['bathy'][...] != 1.)
    
    # volumes.petsc
    vol = ma.array(f['dv'][...], mask=msk)
    print(vol.shape)
    writePETScVector('geometry/volumes.petsc', vol.compressed())




