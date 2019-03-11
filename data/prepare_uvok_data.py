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

class Context():
    '''An empty class to gather and pass script data.'''
    pass

###

def open_mat(file_path):
    '''Open a Matlab file. Use scipy if mat-v5 or h5py if newer.'''
    try:
        f = spio.loadmat(file_path)
        return f
    except:
        f = h5.File(file_path)
        return f

### -----------------------------------------------------------------------------------------------

def writePETScVector(filename, vec):
    '''Write a numpy 1d array into a file in PETSc vector format.'''
    print(filename)
    with open(filename, 'wb') as f:
        # vec id, nvec, data
        nvec, = vec.shape
        np.array([1211214, nvec], dtype='>i4').tofile(f)
        np.array(vec, dtype='>f8').tofile(f)

def writePETScMatrix(filename, mat):
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

def prepare_forcing(ctx):
#    ctx.tmm_path        = 'tmm'
#    ctx.tmm_matlab_path = 'tmm_matlab_code'
#    ctx.uvic_tmm_path   = 'UVic_Kiel_increase_isopyc_diff'
#    ctx.uvic_bgc_path   = 'UVic_Kiel_increase_isopyc_diff_model_data'
    pass
    # prepare boundary
    file_path = ctx.uvic_bgc_path + '/BiogeochemData/UVOK_input_data.mat'
    f = open_mat(file_path)
    for i in range(12):
        var = ma.array(f['aice'][:,:,i].transpose(), mask=ctx.mask[0,:,:])
        var = var.compressed()
        writePETScVector('forcing/boundary/aice_{:02d}.petsc'.format(i), var)



    
    
#    file_path = ctx.uvic_bgc_path + '/grid.mat'
#    f = open_mat(file_path)
## grid.mat
## file ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat
## Hierarchical Data Format (version 5) with 512 bytes user block
#with h5.File('./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat') as f:
#    # forcing/boundary/latitude.petsc
#    y = f['y'][...].transpose()
#    phi = np.zeros((ny, nx))
#    phi[...] = y[:]
#    phi = ma.array(phi, mask=msk[0,:,:])
#    phi = phi.compressed()
#    writePETScVector('forcing/boundary/latitude.petsc', phi)

#    file_path = ctx.uvic_bgc_path + '/BiogeochemData/UVOK_input_data.mat'
#    f = open_mat(file_path)
##    Fe = ma.array(f['Fe'][:,:,:,0], mask=ctx.mask.transpose())
#    Fe = ma.array(f['Fe'][:,:,:,0].transpose(), mask=ctx.mask)
#    print(Fe.shape)
#    print(Fe.flags)
#
##    Fe = Fe.transpose()
#    Fe = Fe.compressed()
#    print(ctx.new_index)
#
#    Fe = Fe[ctx.new_index]
#    writePETScVector('forcing/domain/Fe.petsc', Fe)


#    print(f['aice'][...].shape)
#    print(f['aice'][...].flags)
#    print(f['aice'][...].transpose().shape)
#    print(f['aice'][...].transpose().flags)
#    print(ctx.mask.shape)
#    print(np.ascontiguousarray(f['aice'][...]).shape)
#    aice = f['aice']
#    print(f['aice'][0,:,0])
#    print(aice[0,:,0])
#    aice = ma.zeros(ctx.mask.shape, mask=ctx.mask)
#    aice[0,:,:] = f['aice'][:,:,0].transpose()
#    aice = ma.array(f['aice'][:,:,0].transpose(), mask=ctx.mask[0,:,:])
#    print(aice.shape)
#    print(aice.flags)
#    print(aice[4,:])
#    print(aice[5,:])
#    print(aice[6,:])
#    aice = ma.array(f['aice'][...].transpose()[0,:,:], mask=ctx.mask[0,:,:].transpose())
#    aice = ma.array(f['aice'][:,:,0], mask=ctx.mask[0,:,:].transpose())
#    aice = aice.transpose()
#    aice = aice.compressed()
#    print(aice.shape)
#    print(ctx.new_index.shape)


#    print(aice.shape)
#    print(aice[:100])
#    aice = aice[order]
#    print(aice[:100])


#with h5.File('UVic_Kiel_increase_isopyc_diff/Matrix1/Data/profile_data.mat') as f:
#    Ir_pre = np.array(f['Ir_pre'], dtype='i4') - 1
##    print(f['Ir_pre'][...].shape)
##    print(Ir_pre.shape)
#

#    writePETScVector('forcing/boundary/aice.petsc', aice)

#    import matplotlib
#    matplotlib.use("TkAgg")
#    import matplotlib.pyplot as plt
#
#    plt.imshow(f['aice'][:,:,0])
#    plt.savefig('aice.png')
#
#    plt.imshow(aice[:,:])
#    plt.savefig('aice2.png')

#    plt.spy(ctx.mask[0,:,:].transpose())
#    plt.savefig('mask.png')





#    print(f.keys())
#    print(f['aice'][...].shape)

#    file_path = ctx.uvic_tmm_path + '/Matrix1/Data/boxes.mat'
#    f = open_mat(file_path)
#    print(f.keys())
#    print(f['nb'][...])

#print(spio.whosmat(file_path))
#print(file_path)
#f = open(file_path)
#f.close()
#with h5.File(file_path) as f:
#with h5.File('UVic_Kiel_increase_isopyc_diff_model_data/BiogeochemData/UVOK_input_data.mat') as f:
#    print(f.keys())

#file_path = uvic_bgc_path + '/BiogeochemData/UVOK_input_data.mat'
#uvok_input = spio.loadmat(file_path)
##print(uvok_input.keys())
##dict_keys(['__header__', '__version__', '__globals__', 'Fe', 'aice', 'hice', 'hsno', 'swrad', 'wind'])
#
#for key in uvok_input.keys():
#    print(key, type(uvok_input[key]))
#    if isinstance(uvok_input[key], np.ndarray):
#        print(uvok_input[key].shape)

#load(uvokInputDataFile,'aice')
#load(uvokInputDataFile,'hice')
#load(uvokInputDataFile,'hsno')
#load(uvokInputDataFile,'wind')
#load(uvokInputDataFile,'Fe')
#load(uvokInputDataFile,'swrad')


def prepare_geometry(ctx):
    '''
        `volumes.petsc`
        `landSeaMask.petsc`
    '''
    pass

def prepare_ini(ctx):
    pass

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


def prepare_transport(ctx):
    pass

## Ae_00
## file UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_01.mat
## Hierarchical Data Format (version 5) with 512 bytes user block
#with h5.File('UVic_Kiel_increase_isopyc_diff/Matrix1/TMs/matrix_nocorrection_01.mat') as f:
#    data = f['Aexp']['data'][...]
#    ir = f['Aexp']['ir'][...]
#    jc = f['Aexp']['jc'][...]
#    Ae_ = spsp.csr_matrix((data, ir, jc))
#    Ae_ = Ae_.transpose(copy=True)
##    print(f['Aexp']['jc'][0:10])
##    print(Ae_[330,330])
#    print(Ae_.shape)
#    print(Ae_.nnz)
##    print(Ae_.has_sorted_indices)
#print(Ae_.data[:5])
#print(Ae_[0,:5])
#print(Ae_[1,:5])
#
## Ae index
## file UVic_Kiel_increase_isopyc_diff/Matrix1/Data/profile_data.mat
## Hierarchical Data Format (version 5) with 512 bytes user block
#with h5.File('UVic_Kiel_increase_isopyc_diff/Matrix1/Data/profile_data.mat') as f:
#    Ir_pre = np.array(f['Ir_pre'], dtype='i4') - 1
##    print(f['Ir_pre'][...].shape)
##    print(Ir_pre.shape)
#
#order = Ir_pre[0,:]
#print(order[:10])
#order = np.argsort(order)
#print(order[:10])
#
#Ae = Ae_.tocoo(copy=True)
#Ae = spsp.coo_matrix((Ae.data,(order[Ae.row],order[Ae.col])))
#Ae = Ae.tocsr()
#print(Ae.shape)
#print(Ae.nnz)
##print(Ae.has_sorted_indices)
#print(Ae.data[:5])
#print(Ae[0,:5])
#print(Ae[1,:5])
#
#I = spsp.eye(Ae.shape[0])
#Ae = I + 28800.0*Ae
##Ae = I + (1.0/1095.0)*Ae
##1095/365=3 1/d = 8h = 8*3600s = 28800
#print(Ae.shape)
#print(Ae.nnz)
##print(Ae.has_sorted_indices)
#print(Ae.data[:5])
#print(Ae[0,:5])
#print(Ae[1,:5])
#
#writePETScMatrix('Ae_00.petsc', Ae)



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
    ctx.tmm_path        = 'tmm'
    ctx.tmm_matlab_path = 'tmm_matlab_code'
    ctx.uvic_tmm_path   = 'UVic_Kiel_increase_isopyc_diff'
    ctx.uvic_bgc_path   = 'UVic_Kiel_increase_isopyc_diff_model_data'

    # mask
    file_path = ctx.uvic_bgc_path + '/grid.mat'
    f = open_mat(file_path)
    ctx.mask = (f['bathy'][...] != 1.)
    
    # new_order
    file_path = ctx.uvic_tmm_path + '/Matrix1/Data/profile_data.mat'
    f = open_mat(file_path)
    ctx.new_order = np.array(f['Irr'][0,:] - 1, dtype='i4')

    # new_index
    ctx.new_index = np.array(f['Ir_pre'][0,:] - 1, dtype='i4')

    prepare_forcing(ctx)
    prepare_geometry(ctx)
    prepare_ini(ctx)
    prepare_transport(ctx)

if __name__ == "__main__":
    ctx = Context()
    prepare_uvok_data(ctx)





#    # geometry/landSeaMask.petsc
#    lsm = spsp.csr_matrix(f['ideep'][...])
#    writePETScMatrix('geometry/landSeaMask.petsc', lsm)

#    # geometry/volumes.petsc
#    vol = ma.array(f['dv'][...], mask=msk)
#    vol = np.reshape(vol, (nz, ny*nx)).transpose()
#    vol = np.reshape(vol, (ny, nx, nz))
#    vol = vol.compressed()
#    writePETScVector('geometry/volumes.petsc', vol)

## grid.mat
## file ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat
## Hierarchical Data Format (version 5) with 512 bytes user block
#with h5.File('./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat') as f:
#    # forcing/domain/dz.petsc
#    dz = ma.array(f['dz'][...], mask=msk)
#    dz = np.reshape(dz, (nz, ny*nx)).transpose()
#    dz = np.reshape(dz, (ny, nx, nz))
#    dz = dz.compressed()
#    # scale, m to cm
#    dz = 100*dz
#    writePETScVector('forcing/domain/dz.petsc', dz)












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


#import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.pyplot as plt
#
#n=5000
#plt.spy(Ae_[:n,:n])
#plt.savefig('Ae_.png')
#
#plt.clf()
#plt.spy(Ae[:n,:n])
#plt.savefig('Ae.png')
##plt.show()
#
#Ae_00 = np.fromfile('../../../playground/uvok/metos3d_uvok/data/matrix/Ae_00.petsc', dtype='4>i4,87307>i4,9143571>i4,9143571>f8', count=1)
#
##print(Ae_00.dtype)
##print(Ae_00['f0'].shape)
##print(Ae_00['f1'].shape)
##print(Ae_00['f2'].shape)
##print(Ae_00['f3'].shape)
##print(np.array([0,np.cumsum(Ae_00['f1'])]))
#Ae_00 = spsp.csr_matrix((Ae_00['f3'][0,:],
#                         Ae_00['f2'][0,:],
#                         np.insert(np.cumsum(Ae_00['f1'][0,:]), 0, 0),
#                        ))
#print(Ae_00.shape)
#print(Ae_00.nnz)
##print(Ae.has_sorted_indices)
#print(Ae_00.data[:5])
#print(Ae_00[0,:5])
#print(Ae_00[1,:5])
#
#plt.clf()
#plt.spy(Ae_00[:n,:n])
#plt.savefig('Ae_00.png')

## grid.mat
## file ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat
## Hierarchical Data Format (version 5) with 512 bytes user block
#with h5.File('./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat') as f:
#
##    print(f.keys())
##    print(f['deltaT'][...])
#    ## file ./UVic_Kiel_increase_isopyc_diff_model_data/grid.mat
#    nz, ny, nx = msk.shape
