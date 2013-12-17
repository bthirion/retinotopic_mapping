"""
Script that converts freesurfer outputs (meshes and textures) to .gii files
"""
import commands
import os

#import config_look_loc
import config_retino_7T
#paths = config_look_loc.make_paths(freesurfer=True)
paths = config_retino_7T.make_paths()

# possibly loop on the subjects
subjects = paths.keys()

list_mesh = [ 'rh.pial', 'lh.pial', 'rh.white', 'lh.white', 'rh.inflated',
              'lh.inflated', 'lh.sphere.reg', 'rh.sphere.reg']
list_tex = ['rh.curv', 'lh.curv', 'lh.avg_curv', 'rh.avg_curv']

for subject in subjects:
    # get the common directory of all the stuff
    surf_path = os.path.dirname(paths[subject]['left_mesh'])

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
