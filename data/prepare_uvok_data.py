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
#   This file is based on:
#       https://github.com/samarkhatiwala/tmm/blob/master/models/current/uvok1.0/matlab/make_input_files_for_uvok_model.m
#

import scipy.sparse as spsp
import scipy.io as spio
import h5py as h5
import numpy as np
import numpy.ma as ma

### -----------------------------------------------------------------------------------------------

class Context():
    '''An empty class to gather and pass script data.'''
    pass

### -----------------------------------------------------------------------------------------------

def open_mat(file_path):
    '''Open a matlab file. Use scipy if mat-v5 or h5py if newer.'''
    print('Open Matlab file ...', file_path)
    try:
        f = spio.loadmat(file_path)
        return f
    except:
        f = h5.File(file_path)
        return f

### -----------------------------------------------------------------------------------------------

def write_petsc_vector_list(filename, vec_list):
    '''Write a list of vectors.'''
    for i in range(len(vec_list)):
        filename_i = filename + '_{:02d}.petsc'.format(i)
        write_petsc_vector(filename_i, vec_list[i])

def write_petsc_vector(filename, vec):
    '''Write a numpy 1d array into a file in PETSc vector format.'''
    print(filename)
    with open(filename, 'wb') as f:
        # vec id, nvec, data
        nvec, = vec.shape
        np.array([1211214, nvec], dtype='>i4').tofile(f)
        np.array(vec, dtype='>f8').tofile(f)

def write_petsc_matrix(filename, mat):
    '''Write a scipy sparse matrix into a file in PETSc sparse format.'''
    print(filename)
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

def read_from_mat_file(ctx, filename, varname):
    '''Read an (arbitrarily shaped) variable from a matlab file.'''
    f = open_mat(filename)
    var = f[varname][...]
    return var

def read_2d_list_from_mat_file(ctx, filename, varname, n):
    '''Read an array of 2d-arrays from a matlab file. Return as a list of arrays.
        Apply a grid mask to and compress each 2d-array.'''
    f = open_mat(filename)
    var_list = []
    for i in range(n):
        var = ma.array(f[varname][:,:,i].transpose(), mask=ctx.mask[0,:,:])
        var = var.compressed()
        var_list.append(var)
    return var_list

def read_3d_list_from_mat_file(ctx, filename, varname, n):
    '''Read an array of 3d-arrays from a matlab file. Return as a list of arrays.
        Apply a grid mask to and compress each 3d-array.
        Re-order the compressed array.'''
    f = open_mat(filename)
    var_list = []
    for i in range(n):
        var = ma.array(f[varname][:,:,:,i].transpose(), mask=ctx.mask)
        var = var.compressed()
        var = var[ctx.new_index]
        var_list.append(var)
    return var_list

### -----------------------------------------------------------------------------------------------

def reshape_and_compress(ctx, var):
    '''Create a masked array from variable. Reshape from leading x to leading z dimension.
        Preserve xy order. Compress.'''
    var = ma.array(var, mask=ctx.mask)
    nz, ny, nx = ctx.mask.shape
    var = np.reshape(var, (nz, ny*nx)).transpose()
    var = np.reshape(var, (ny, nx, nz))
    var = var.compressed()
    return var

def reorder_matrix(ctx, mat):
    '''Create a csr sparse matrix.'''
    mat = spsp.csr_matrix((mat['data'], mat['ir'], mat['jc']))
    mat = mat.transpose()
    mat = mat.tocoo()
    mat = spsp.coo_matrix((mat.data,(ctx.new_order[mat.row],ctx.new_order[mat.col])))
    mat = mat.tocsr()
    return mat

### -----------------------------------------------------------------------------------------------

def prepare_forcing_boundary(ctx, file_in, var_name, n, file_name):
    '''Prepare boundary (2d) data from given file.'''
    var = read_2d_list_from_mat_file(ctx, file_in, var_name, n)
    file_out = 'forcing/boundary/' + file_name
    write_petsc_vector_list(file_out, var)

def prepare_forcing_domain(ctx, file_in, var_name, n, file_name):
    '''Prepare domain (3d) data from given file.'''
    var = read_3d_list_from_mat_file(ctx, file_in, var_name, n)
    file_out = 'forcing/domain/' + file_name
    write_petsc_vector_list(file_out, var)

def prepare_forcing(ctx):
    '''Prepare boundary and domain forcing.'''
    # prepare boundary
    # latitude
    file_in = ctx.uvic_tmm_path + '/Matrix1/Data/boxes.mat'
    lat = read_from_mat_file(ctx, file_in, 'Ybox')
    lat = lat[0,:6386]
    write_petsc_vector('forcing/boundary/latitude.petsc', lat)

    file_in = ctx.uvic_bgc_path + '/BiogeochemData/UVOK_input_data.mat'
    prepare_forcing_boundary(ctx, file_in, 'aice', 12, 'aice')
    prepare_forcing_boundary(ctx, file_in, 'hice', 12, 'hice')
    prepare_forcing_boundary(ctx, file_in, 'hsno', 12, 'hsno')
    prepare_forcing_boundary(ctx, file_in, 'wind', 12, 'wind')
    prepare_forcing_boundary(ctx, file_in, 'swrad', 12, 'swrad')

    # emp
    file_in = ctx.uvic_bgc_path + '/grid.mat'
    da = read_from_mat_file(ctx, file_in, 'da')
    da = ma.array(da, mask=ctx.mask)
    da = da[0,:,:].compressed()
    #
    file_in = ctx.uvic_bgc_path + '/GCM/FreshWaterForcing_gcm.mat'
    emp = read_2d_list_from_mat_file(ctx, file_in, 'EmPgcm', 12)
    emp = np.array(emp)
    emp = emp - np.mean(emp.dot(da))/np.sum(da)
    emp = 100*emp
    emp = [emp[i,:] for i in range(emp.shape[0])]
    #
    file_out = 'forcing/boundary/' + 'EmP'
    write_petsc_vector_list(file_out, emp)
    
    # prepare domain
    file_in = ctx.uvic_bgc_path + '/BiogeochemData/UVOK_input_data.mat'
    prepare_forcing_domain(ctx, file_in, 'Fe', 12, 'Fe_dissolved')
    # salt
    file_in = ctx.uvic_bgc_path + '/GCM/Salt_gcm.mat'
    salt_list = read_3d_list_from_mat_file(ctx, file_in, 'Sgcm', 12)
    salt_list = [(s - 35.0)/1000.0 for s in salt_list]
    write_petsc_vector_list('forcing/domain/Ss', salt_list)
    # temp
    file_in = ctx.uvic_bgc_path + '/GCM/Theta_gcm.mat'
    prepare_forcing_domain(ctx, file_in, 'Tgcm', 12, 'Ts')
    # dz
    file_in = ctx.uvic_bgc_path + '/grid.mat'
    dz = read_from_mat_file(ctx, file_in, 'dz')
    dz = reshape_and_compress(ctx, dz)
    dz = 100*dz
    write_petsc_vector('forcing/domain/dz.petsc', dz)
    # zt
    file_in = ctx.uvic_tmm_path + '/Matrix1/Data/boxes.mat'
    zt = read_from_mat_file(ctx, file_in, 'Zbox')
    zt = zt[0,ctx.new_index]
    zt = 100*zt
    write_petsc_vector('forcing/domain/zt.petsc', zt)

def prepare_geometry(ctx):
    '''Prepare geometry files.'''
    # land sea mask
    file_in = ctx.uvic_bgc_path + '/grid.mat'
    lsm = read_from_mat_file(ctx, file_in, 'ideep')
    lsm = spsp.csr_matrix(lsm)
    write_petsc_matrix('geometry/landSeaMask.petsc', lsm)
    # volumes
    vol = read_from_mat_file(ctx, file_in, 'dv')
    vol = reshape_and_compress(ctx, vol)
    write_petsc_vector('geometry/volumes.petsc', vol)

def prepare_ini(ctx):
    '''Prepare initial tracer concentrations.'''
    names = ['dic','c14','alk','o2','po4','phyt','zoop','detr','no3','diaz']
    for name in names:
        file_in = ctx.uvok_path + '/matlab/InitialConditionProfiles/' + name + '.dat'
        vec1d = np.loadtxt(file_in)[:,1]
        vec3d = np.zeros(ctx.mask.shape)
        vec3d[...] = vec1d[:, np.newaxis, np.newaxis]
        vec3d = reshape_and_compress(ctx, vec3d)
        # write
        file_out = 'ini/' + name + 'ini.petsc'
        write_petsc_vector(file_out, vec3d)

def prepare_transport(ctx):
    '''Prepare transport matrices.'''
    for i in range(12):
        file_in = ctx.uvic_tmm_path + '/Matrix1/TMs/matrix_nocorrection_{:02d}.mat'.format(i+1)
        f = open_mat(file_in)
        # exp
        Ae = reorder_matrix(ctx, f['Aexp'])
        I = spsp.eye(Ae.shape[0])
        Ae = I + 28800.0*Ae
        file_out = 'transport/Ae_{:02d}.petsc'.format(i)
        write_petsc_matrix(file_out, Ae)
        # imp
        Ai = reorder_matrix(ctx, f['Aimp'])
        file_out = 'transport/Ai_{:02d}.petsc'.format(i)
        write_petsc_matrix(file_out, Ai)

### -----------------------------------------------------------------------------------------------

def prepare_uvok_data(ctx):
    '''
        This routine prepares the required data for UVOK simulation.
        Assume we are located in data/ and have the following directory structure:
        
        data/
        ├── forcing
        │   ├── boundary
        │   └── domain
        ├── geometry
        ├── ini
        └── transport

        We also prepared the TMM/UVIC/UVOK sources as described in `../README.md`.
        '''
    # input directories
    ctx.uvok_path       = 'tmm/models/current/uvok1.0'
    ctx.uvic_tmm_path   = 'UVic_Kiel_increase_isopyc_diff'
    ctx.uvic_bgc_path   = 'UVic_Kiel_increase_isopyc_diff_model_data'
    # mask
    file_path = ctx.uvic_bgc_path + '/grid.mat'
    f = open_mat(file_path)
    ctx.mask = (f['bathy'][...] != 1.)
    # new_order, new_index
    file_path = ctx.uvic_tmm_path + '/Matrix1/Data/profile_data.mat'
    f = open_mat(file_path)
    ctx.new_order = np.array(f['Irr'][0,:] - 1, dtype='i4')
    ctx.new_index = np.array(f['Ir_pre'][0,:] - 1, dtype='i4')

    prepare_forcing(ctx)
    prepare_geometry(ctx)
    prepare_ini(ctx)
    prepare_transport(ctx)

### -----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    ctx = Context()
    prepare_uvok_data(ctx)

