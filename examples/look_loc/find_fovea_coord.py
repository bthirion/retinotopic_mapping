"""
This scipt aims at finding the coordinates fo the fovea in the subjects
of this dataset, so that it can be sued in future experiments.

Note
----
First run export SUBJECTS_DIR=/volatile/thirion/data/lookloc

Author: Bertrand Thirion, 2012
"""
import os.path as op
import numpy as np
from nibabel.gifti import read
try:
    from parietal.surface_operations.mesh_processing import mesh_arrays
except:
    from retino.mesh_processing import mesh_arrays
try:
    from mayavi import mlab
except:
    from enthought.mayavi import mlab

from retino.angular_analysis import load_texture
from retino.group_analysis import resample_to_average
from config_look_loc import make_paths

# get all the ,necessary paths
paths = make_paths()

# take the fsaverage model
left_ref_mesh = op.join(paths[paths.keys()[0]]['base'],
                        '../fsaverage/lh.white.gii')
left_vertices, _ = mesh_arrays(left_ref_mesh)
right_ref_mesh = op.join(paths[paths.keys()[0]]['base'],
                        '../fsaverage/rh.white.gii')
right_vertices, _ = mesh_arrays(right_ref_mesh)


### fixme: this is copied from group_analysis
left_pos, right_pos = [], []
for subject in paths.keys():
    print subject
    func_path = op.join(paths[subject]['base'], paths[subject]['acquisition'],
                        'analysis')
    # set all the paths
    ltex_path = op.join(func_path, 'left_phase_ring.gii')
    rtex_path = op.join(func_path, 'right_phase_ring.gii')

    for side in ['left', 'right']:
        #resample on tyhe average
        tex_path = ltex_path if side == 'left' else rtex_path
        resampled = resample_to_average(tex_path, paths[subject]['fs_subj'],
                                        side)
        # mask the values that are small enough
        rring = load_texture(resampled).ravel()
        mask = rring < rring.min() + 1.
        if side == 'left':
            left_pos.append(left_vertices[mask].mean(0))
        else:
            right_pos.append(right_vertices[mask].mean(0))

left_pos = np.array(left_pos)
right_pos = np.array(right_pos)
print left_pos.std(0), right_pos.std(0)

if True:
    left_fovea_, right_fovea_ = np.median(left_pos, 0), np.median(right_pos, 0)
else:
    left_fovea_, right_fovea_ = np.mean(left_pos, 0), np.mean(right_pos, 0)

rescaling = 1.03
# find the nearest node on the template mesh
left_fovea = left_vertices[np.argmin(
        np.sum((left_vertices - left_fovea_ * rescaling) ** 2, 1))]

right_fovea = right_vertices[np.argmin(
        np.sum((right_vertices - right_fovea_ * rescaling) ** 2, 1))]

# Plotting to look at the result

mlab.figure(bgcolor=(.8, .8, .8))
for point in [left_fovea, right_fovea]:
    x, y, z = point[:, np.newaxis]
    pts = mlab.points3d(x, y, z, 5,
                        scale_factor=5., resolution=10, scale_mode='none')

for mesh in [left_ref_mesh, right_ref_mesh]:
    vertices = read(mesh).darrays[0].data
    xv, yv, zv = vertices.T
    triangles = read(mesh).darrays[1].data
    mlab.triangular_mesh(xv, yv, zv, triangles, transparent=False,
                         opacity=1., name='', color=(0.5, 0.5, 0.5))

mlab.show()
