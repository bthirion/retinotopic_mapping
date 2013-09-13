"""
Script that converts freesurfer outputs (meshes and textures) to .gii files
This is necessary for the sake of vizualization with standard viewers.

Author: Bertrand Thirion, 2012
"""
import commands
import os
#import config_retino_7T

#_, data_path, subject_info = config_retino_7T.init_config()
 
#subjects = subject_info.keys()

# data directory
main_dir = '/neurospin/tmp/retino/7T/'

list_mesh = [ 'rh.pial', 'lh.pial', 'rh.white', 'lh.white', 'rh.inflated',
              'lh.inflated', 'lh.sphere.reg', 'rh.sphere.reg']
list_tex = ['rh.curv', 'lh.curv', 'lh.avg_curv', 'rh.avg_curv']

for subject in ['dummy_1mm']:#subjects:
    if subject == 'gm110134': 
        continue
    # get the common directory of all the stuff
    surf_path = os.path.join(main_dir, subject, 't1', subject, 'surf')

    # convert the meshes
    for fname in list_mesh:
        input_path = os.path.join(surf_path, fname)
        output_path = input_path + '.gii'
        cmd = 'mris_convert %s %s' % (input_path, output_path)
        print commands.getoutput(cmd)

    # convert the textures
    for fname in list_tex:
        input_path = os.path.join(surf_path, fname)
        output_path = input_path + '.gii'
        cmd = 'mri_convert %s %s' % (input_path, output_path)
        print commands.getoutput(cmd)
