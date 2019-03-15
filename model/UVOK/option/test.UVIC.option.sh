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

# debug
-Metos3DDebugLevel                                  3

# geometry
-Metos3DGeometryType                                Profile
-Metos3DProfileInputDirectory                       data/geometry/
-Metos3DProfileMaskFile                             landSeaMask.petsc
-Metos3DProfileVolumeFile                           volumes.petsc

# bgc tracer
-Metos3DTracerCount                                 10
-Metos3DTracerName                                  DIC,C14,ALK,O2,PO4,PHY,ZOO,DET,NO3,DIAZ
-Metos3DTracerInputDirectory                        data/ini/
-Metos3DTracerInitFile                              dicini.petsc,c14ini.petsc,alkini.petsc,o2ini.petsc,po4ini.petsc,phytini.petsc,zoopini.petsc,detrini.petsc,no3ini.petsc,diazini.petsc
-Metos3DTracerOutputDirectory                       work/
-Metos3DTracerOutputFile                            dic.petsc,c14.petsc,alk.petsc,o2.petsc,po4.petsc,phyt.petsc,zoop.petsc,detr.petsc,no3.petsc,diaz.petsc
# weight with volumes and sum up
-Metos3DTracerMonitor

# diagnostic variables
-Metos3DDiagnosticCount                             0
# weight with volumes and sum up
#-Metos3DDiagnosticMonitor

# bgc parameter
-Metos3DParameterCount                              0
#-Metos3DParameterValue                              

# bgc boudary conditions
-Metos3DBoundaryConditionCount                      7
-Metos3DBoundaryConditionInputDirectory             data/forcing/boundary/
-Metos3DBoundaryConditionName                       Latitude,Wind,SWRAD,IceArea,IceHeight,SnowHeight,EmP
# latitude
# wind
# short wave radiation
# ice area
# ice height
# snow height
# evaporation minus precipitation
-Metos3DLatitudeCount                               1
-Metos3DLatitudeFileFormat                          latitude.petsc
-Metos3DWindCount                                   12
-Metos3DWindFileFormat                              wind_$02d.petsc
-Metos3DSWRADCount                                  12
-Metos3DSWRADFileFormat                             swrad_$02d.petsc
-Metos3DIceAreaCount                                12
-Metos3DIceAreaFileFormat                           aice_$02d.petsc
-Metos3DIceHeightCount                              12
-Metos3DIceHeightFileFormat                         hice_$02d.petsc
-Metos3DSnowHeightCount                             12
-Metos3DSnowHeightFileFormat                        hsno_$02d.petsc
-Metos3DEmPCount                                    12
-Metos3DEmPFileFormat                               EmP_$02d.petsc

# bgc domain conditions
-Metos3DDomainConditionCount                        5
-Metos3DDomainConditionInputDirectory               data/forcing/domain/
-Metos3DDomainConditionName                         WithinLayer,LayerHeight,Temperature,Salinity,FeDiss
# points within layer
# layer height
# temperature
# salinity
# dissolved iron
-Metos3DWithinLayerCount                            1
-Metos3DWithinLayerFileFormat                       zt.petsc
-Metos3DLayerHeightCount                            1
-Metos3DLayerHeightFileFormat                       dz.petsc
-Metos3DTemperatureCount                            12
-Metos3DTemperatureFileFormat                       Ts_$02d.petsc
-Metos3DSalinityCount                               12
-Metos3DSalinityFileFormat                          Ss_$02d.petsc
-Metos3DFeDissCount                                 12
-Metos3DFeDissFileFormat                            Fe_dissolved_$02d.petsc

# transport
-Metos3DTransportType                               Matrix
-Metos3DMatrixInputDirectory                        data/transport/
#-Metos3DMatrixCount                                 12
-Metos3DMatrixCount                                 2
-Metos3DMatrixExplicitFileFormat                    Ae_$02d.petsc
-Metos3DMatrixImplicitFileFormat                    Ai_$02d.petsc

# time stepping
-Metos3DTimeStepStart                               0.0
#-Metos3DTimeStepCount                               1095
-Metos3DTimeStepCount                               1
-Metos3DTimeStep                                    0.0009132420091324

# solver
-Metos3DSolverType                                  Spinup
#-Metos3DSpinupCount                                 3000
-Metos3DSpinupCount                                 1
-Metos3DSpinupMonitor

## solver (PETSc)
#-Metos3DSolverType                                  Newton
#-Metos3DNewton_snes_type                            ls
#-Metos3DNewton_snes_view
#-Metos3DNewton_snes_ksp_ew
#-Metos3DNewton_snes_monitor
#-Metos3DNewton_snes_linesearch_monitor
#-Metos3DNewton_ksp_type                             gmres
#-Metos3DNewton_ksp_monitor
#-Metos3DNewton_ksp_view


