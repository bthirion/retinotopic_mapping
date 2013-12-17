"""
Fit the visual areas to individual subjects data

Note
====
first type export SUBJECTS_DIR=/volatile/thirion/data/lookloc/

Todo
====

Author:Bertrand Thirion, 2012
"""

import numpy as np
import os.path as op
import enthought.mayavi.mlab as mlab
from parietal.glm_files_layout.cortical_glm import load_texture
from config_look_loc import make_paths
from angular_analysis import (phase_maps, delineate_areas)
from visu_mlab import plot_retino_image
try:
    from parietal.surface_operations.mesh_processing import isomap_patch
except:
    from parietal_copy.mesh_processing import isomap_patch


#########################################################################
# Main script
#########################################################################

paths = make_paths()
subjects = paths.keys()
all_reg = ['sin_ring_pos', 'cos_ring_pos', 'sin_ring_neg', 'cos_ring_neg',
           'sin_wedge_neg', 'cos_wedge_neg', 'sin_wedge_pos', 'cos_wedge_pos']
offset_ring, offset_wedge = np.pi, - np.pi / 2
do_phase_unwrapping = True
do_plot = False

for subject in subjects:
    print subject
    tex_dir = op.join(paths[subject]['base'], 
                      paths[subject]['acquisition'], 'analysis')

    for side in ['right', 'left']:
        # First compute a mask of the active regions, by thresholding the
        # so-called "average stat" map
        mesh = paths[subject]['%s_inflated' % side]
        mask = load_texture(op.join(tex_dir, side + '_mask.gii')) > 0

        # Then, given the mask, extract the beta data from each subject
        data = dict([ (r, load_texture(
                        op.join(tex_dir, '%s_%s_con.gii' % (side, r))).ravel()
                       [mask]) for r in all_reg])
    
        # Compute the phase value of the data
        phase_wedge, phase_ring, hemo = phase_maps(
        data, offset_ring, offset_wedge, do_ring=True, do_wedge=True, 
        do_phase_unwrapping=do_phase_unwrapping, mesh=mesh, mask=mask)

        # get planar coordinates on the surface
        planar_coord = isomap_patch(mesh, mask)

        # template-based visual area delineation
        visual_areas = delineate_areas(phase_ring, phase_wedge,  hemo, mesh, 
                                       mask, planar_coord, side, tex_dir)
    
        #--------------------------------------------------------------------
        # plot the result
        #--------------------------------------------------------------------
        if do_plot:
            wedge_tex = op.join(tex_dir, side + '_phase_wedge.gii')
            curv_tex = op.join(op.dirname(mesh), '%sh.avg_curv.gii' % side[0])
            
            f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
            mlab.clf()
            if side == 'left':
                vmin, vmax = -3.14, .5
            else:
                vmin, vmax = -.5, np.pi
            plot_retino_image(mesh, name="%s_hemisphere" %side, tf=None, 
                              curv=load_texture(curv_tex), mask=mask,
                              tex=load_texture(wedge_tex), vmin=vmin, vmax=vmax)
            mlab.view(280, 120)
            f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
            mlab.clf()
            plot_retino_image(mesh, name="%sHemisphere" % side, tf=None, 
                              mask=mask, tex=visual_areas, 
                              curv=load_texture(curv_tex), vmin=0, vmax=6)
            mlab.view(280, 120)
    
