"""
Visualization of the results of the retinotopy.

Author: Bertrand Thirion, 2012

"""

import numpy as np
import os.path as op

import enthought.mayavi.mlab as mlab
from nibabel.gifti import read

from retino.visu_mlab import plot_retino_image
import config_retino_7T

# get the paths information
_, main_dir, subject_info = config_retino_7T.init_config()

sides = ['left', 'right']
for subject in ['eb120536']:#subject_info.keys():
    if subject == 'gm110134':
        continue # this one does not work
    print subject

    # set all the paths
    func_path = op.join(main_dir, subject, 'fmri/analysis')
    anat_path = op.join(main_dir, subject, 't1', subject, 'surf')
    for side in sides:
        mesh_path = op.join(anat_path, '%sh.white.gii' % side[0])
        tex_path = op.join(func_path, '%s_phase_wedge.gii' % side)
        inflated_mesh_path = op.join(anat_path, '%sh.inflated.gii' % side[0])
        curv_path = op.join(anat_path, '%sh.avg_curv.gii' % side[0])
    
        # get the masks
        mask_path = op.join(func_path, '%s_mask.gii' % side)

        # read the texture and curvature data
        tex = read(tex_path).darrays[0].data.ravel()
        curv = (read(curv_path).darrays[0].data < 0).astype(np.int)
        mask = read(mask_path).darrays[0].data

        ### Plot meshes
        f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
        mlab.clf()
        vmin, vmax = -np.pi, np.pi
        plot_retino_image(inflated_mesh_path, 
                          name="%s hemisphere" % side, 
                          mask=mask,
                          tf=None, 
                          tex=tex, 
                          curv=curv, 
                          vmin=vmin, vmax=vmax)
        if side == 'left':
            mlab.view(280, 120)
        else:
            mlab.view(250, 120)
        mlab.savefig(op.join(func_path, '%s_%s.png') % (subject, side))

mlab.show()
