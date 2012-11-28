"""
Several helpers for visualization of brain meshes and activity maps
for retinotopic studies.

Author: Virgile Fritsch, Bertrand Thirion 2010-2011
"""

import numpy as np

import enthought.mayavi.mlab as mlab

try:
    from parietal.surface_operations.mesh_processing import mesh_arrays
except:
    from parietal_copy.mesh_processing import mesh_arrays


def plot_retino_image(mesh_path, name, tf=None, tex=None, curv=None, mask=None,
                     vmin=0, vmax=0):
    """Custom function to plot mesh in the case of retinotopic mapping  

    Parameters
    ----------
    mesh_path: string, path of the mesh to plot
    name: string, object identifier
    tf: (4, 4)  array
        affine transformation that has to be applied to the mesh
    tex: array of shape (mes.n_vertices),
         texture to plot on the surface
    curv: array of shape (mes.n_vertices),
          curvature texture to plot on the surface
    mask: array of shape (mes.n_vertices),
          A mask for functional regions of interest
    vmin, vmax: bounds for color clipping
    """
    vertices, triangles =  mesh_arrays(mesh_path)
    af_coord = np.hstack((vertices, np.ones((vertices.shape[0], 1))))
    if tf is not None:
        af_coord = np.dot(af_coord, tf.T)
    
    vertices = af_coord[:, :3]
    x, y, z = vertices.T

    # it is sexpected that the curvature takes values between 0 and 1
    if curv is not None:
        cmin = 2 * curv.min() - curv.max()
        cmax = 2 * curv.max() - curv.min()
        mlab.triangular_mesh(x, y, z, triangles, transparent=False, opacity=1.,
                             name=name, scalars=curv, colormap="bone",
                             vmin=cmin, vmax=cmax)
    else:
        mlab.triangular_mesh(x, y, z, triangles, transparent=False, opacity=1.,
                             name=name)
    if tex is not None:
        if mask is not None:
            tex[mask == 0] = vmin - 1
            print np.sum(tex < vmin)
        func_mesh = mlab.pipeline.triangular_mesh_source(x, y, z,
                                                  triangles,
                                                  scalars=tex)
        thresh = mlab.pipeline.threshold(func_mesh, low=vmin)
        mlab.pipeline.surface(thresh, colormap="jet", vmin=vmin,
                              vmax=vmax, transparent=True,
                              opacity=.8)
        mlab.scalarbar(thresh)


def plot_mesh(mesh, name, tf=None, tex=None):
    """ old version -- probably deprecated
    """
    vertices, triangles =  mesh_arrays(mesh)
    af_coord = np.hstack((vertices, np.ones((vertices.shape[0], 1))))
    if tf is not None:
        af_coord = np.dot(af_coord, tf.T)
    
    vertices = af_coord[:, :3]
    x, y, z = vertices.T
    # show with mayavi

    if tex is None:
        mlab.triangular_mesh(x, y, z, triangles, transparent=False, opacity=1.,
                             name=name, color=(0.5, 0.5, 0.))
    else:
        mlab.triangular_mesh(x, y, z, triangles, transparent=False, opacity=1.,
                             name=name, scalars=tex)
