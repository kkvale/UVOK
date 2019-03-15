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
                # link
                # TODO check if present, remove first before link
                try:
                    os.symlink(file_path, os.path.basename(file_path))
                except:
                    pass
                break


