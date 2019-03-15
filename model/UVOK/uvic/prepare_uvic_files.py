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
#       https://github.com/samarkhatiwala/tmm/blob/master/models/current/uvok1.0/src/Makefile
#

import os

base_dir = '.'

search_dirs = [
               'Kiel_Jan_2017_updates_to_UVIC2.9/source',
               'UVic2.9/updates/02/source',
               'UVic2.9/source',
               ]

search_sub_dirs = [
                   'common',
                   'mom',
                   'embm',
                   ]

dir_list = [sdir + '/' + ssdir  for sdir in search_dirs for ssdir in search_sub_dirs]

source_file_list = [
                    'gasbc.F',
                    'gosbc.F',
                    'setvbc.F',
                    'tracer.F',
                    'co2calc.F',
                    'npzd_src.F',
                    'UVic_ESCM.F',
                    ]

header_file_list = [
                    # uvok_calc.F
                    'size.h',
                    'param.h',
                    'pconst.h',
                    'stdunits.h',
                    'coord.h',
                    'csbc.h',
                    'grdvar.h',
                    'levind.h',
                    'mw.h',
                    'scalar.h',
                    'tmngr.h',
                    'npzd.h',
                    'calendar.h',
                    'diaga.h',
                    'ice.h',
                    'atm.h',
                    # gasbc.F
                    'switch.h',
                    'cembm.h',
                    'insolation.h',
                    'solve.h',
                    # tracer.F
                    'accel.h',
                    'cregin.h',
                    'emode.h',
                    'hmixc.h',
                    'timeavgs.h',
                    'vmixc.h',
                    'isopyc.h',
                    'fdift.h',
                    # UVic_ESCM.F
                    'cnep.h',
                    'cprnts.h',
                    'diag.h',
                    'fwa.h',
                    'iounit.h',
                    'mtlm.h',
                    'sed.h',
                    'stab.h',
                    'veg.h',
                    ]

fmt = '{:.<24} '

if __name__ == '__main__':
    
    # source
    print('source:')
    for file in source_file_list:
        print(fmt.format(file), end='')
        for sdir in dir_list:
            file_path = base_dir + '/' + sdir + '/' + file
#            print(file_path, os.path.exists(file_path))
            if os.path.exists(file_path):
                print(file_path)
                # prepend include statement
                cmd = 'cat include_uvok_tmm_options.h {} > {}'.format(file_path, os.path.basename(file_path))
                os.system(cmd)
                break

    # header
    print('header:')
    for file in header_file_list:
        print(fmt.format(file), end='')
        for sdir in dir_list:
            file_path = base_dir + '/' + sdir + '/' + file
#            print(file_path, os.path.exists(file_path))
            if os.path.exists(file_path):
                print(file_path)
                # copy
                cmd = 'cp {} {}'.format(file_path, os.path.basename(file_path))
                os.system(cmd)
                break


