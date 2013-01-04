"""
Visualization of retinotopy results

Author: Bertrand Thirion, 2012
"""

import numpy as np
import os.path as op

import enthought.mayavi.mlab as mlab
from nibabel.gifti import read

from retino.visu_mlab import plot_retino_image
from config_look_loc import make_paths
paths = make_paths()

for subject in paths.keys():
    print subject
    func_path = op.join(paths[subject]['base'], paths[subject]['acquisition'],
                        'analysis')
    # set all the paths
    ltex_path = op.join(func_path, 'left_phase_wedge.gii')
    rtex_path = op.join(func_path, 'right_phase_wedge.gii')
    
    lmesh_path_inflated = paths[subject]['left_inflated']
    lcurv_path = op.join(op.dirname(lmesh_path_inflated), 'lh.avg_curv.gii')
    rmesh_path_inflated = paths[subject]['right_inflated']
    rcurv_path = op.join(op.dirname(lmesh_path_inflated), 'rh.avg_curv.gii')
    
    lmask_path = op.join(func_path, 'left_mask.gii')
    rmask_path = op.join(func_path, 'right_mask.gii')

    # read the texture and curvature data
    ltex = read(ltex_path).darrays[0].data.ravel()
    lcurv = (read(lcurv_path).darrays[0].data < 0).astype(np.int)
    lmask = read(lmask_path).darrays[0].data
    rtex = read(rtex_path).darrays[0].data.ravel()
    rcurv = (read(rcurv_path).darrays[0].data < 0).astype(np.int)
    rmask = read(rmask_path).darrays[0].data

    ### Plot meshes
    
    # left hemisphere
    f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
    mlab.clf()
    plot_retino_image(lmesh_path_inflated, name="LeftHemisphere", mask=lmask,
                      tf=None, tex=ltex, curv=lcurv, vmin=-np.pi, vmax=np.pi)
    #plot_retino_image(lmesh_path_inflated, name="LeftHemisphere", mask=lmask,
    #                  tf=None, tex=ltex, curv=lcurv, vmin=-np.pi, vmax=0)
    mlab.view(280, 120)
    mlab.savefig(op.join(func_path,'%s_%s.png') % (subject, 'left'))

    # right hemisphere
    f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
    mlab.clf()
    plot_retino_image(rmesh_path_inflated, name="RightHemisphere", mask=rmask,
                      tf=None, tex=rtex, curv=rcurv, vmin=-np.pi, vmax=np.pi)
    #plot_retino_image(rmesh_path_inflated, name="RightHemisphere", mask=rmask,
    #                  tf=None, tex=rtex, curv=rcurv, vmin=0, vmax=np.pi)
    mlab.view(250, 120)
    mlab.savefig(op.join(func_path, '%s_%s.png') % (subject, 'right'))
    mlab.show()

    
