"""
Script to perform angular operations on brain activation maps.
This assumes that the standard file layout defined in nipy has been enforced.

Author: Bertrand Thirion, 2010-2011
"""

import numpy as np
import os.path as op
from nibabel import load, save, Nifti1Image
from nipy.labs.spatial_models.discrete_domain import \
     grid_domain_from_binary_array
from nipy.algorithms.graph import wgraph_from_coo_matrix

try:
    from parietal.surface_operations import mesh_processing as mep
except:
    import mesh_processing as  mep

NEGINF = -np.inf

##################################################################
# Ancillary functions
##################################################################

def load_texture(path):
    """Return an array from texture data stored in a gifti file

    Parameters
    ----------
    path string or list of strings
         path of the texture files

    Returns
    -------
    data array of shape (nnode) or (nnode, len(path)) 
         the corresponding data
    """
    from nibabel.gifti import read

    # use alternative libraries than nibabel if necessary
    if hasattr(path, '__iter__'):
        tex_data = []
        for f in path:
            ftex = read(f).getArraysFromIntent('NIFTI_INTENT_TIME_SERIES')
            tex = np.array([f.data for f in ftex])
            tex_data.append(tex)
        tex_data = np.array(tex_data)
        if len(tex_data.shape) > 2:
            tex_data = np.squeeze(tex_data)
    else:
        tex_ = read(path)
        tex_data = np.array([darray.data for darray in tex_.darrays])
    return tex_data


def save_texture(path, data, intent='none', verbose=False):
    """
    volume saving utility for textures
    
    Parameters
    ----------
    path, string, output image path
    data, array of shape (nnode)
          data to be put in the volume
    intent: string, optional
            intent

    Fixme
    -----
    Missing checks
    Handle the case where data is multi-dimensional ? 
    """
    from nibabel.gifti import write, GiftiDataArray, GiftiImage
    if verbose:
        print 'Warning: assuming a float32 gifti file'
    darray = GiftiDataArray().from_array(data.astype(np.float32), intent)
    img = GiftiImage(darrays=[darray])
    write(img, path)


def cc_mask(graph, mask, size_threshold):
    """Filter the mask to keep only connected components above a given size"""
    mg = graph.subgraph(mask.ravel())
    u = mg.cc()
    for k in np.unique(u):
        if np.sum(u == k) < size_threshold:
            u[u == k] = -1
    mask[mask] = u > - 1
    return mask


def cc_mesh_mask(mesh, mask, size_threshold):
    """Filter the mask to keep only connected components above a given size"""
    return cc_mask(mep.mesh_to_graph(mesh), mask, size_threshold)


def cc_array_mask(mask, size_threshold):
    """Filter the mask to keep only connected components above a given size"""
    graph = wgraph_from_coo_matrix(grid_domain_from_binary_array(mask).topology)
    mask[mask] = cc_mask(graph, mask[mask], size_threshold)
    return mask


def main_cc(mesh, mask):
    """Filter the mask to keep only the main connected component"""
    mg = mep.mesh_to_graph(mesh).subgraph(mask.ravel())
    u = mg.cc()
    size_u = np.array([np.sum(u == k) for k in np.unique(u) if k > -1])
    mask[mask] = (u == size_u.argmax())
    return mask


def combine_phase(phase_pos, phase_neg, offset=0, hemo=None):
        """ Combine the phases estimated in two directions"""
        if hemo == None:
            # estimate hemodynamic delay
            hemo = 0.5 * (phase_pos + phase_neg)
            hemo += np.pi * (hemo < - np.pi / 4)
            hemo += np.pi * (hemo < - np.pi / 4)
        
        # first phase estimate
        pr1 = phase_pos - hemo
        pr2 = hemo - phase_neg
        pr2[pr1 - pr2 > np.pi] += 2 * np.pi
        pr1[pr2 - pr1 > np.pi] += 2 * np.pi
        phase = 0.5 * (pr1 + pr2)
        
        # add the offset and bring back to [-pi, +pi]
        phase += offset
        phase += 2 * np.pi * (phase < - np.pi)
        phase -= 2 * np.pi * (phase > np.pi)
        phase += 2 * np.pi * (phase < - np.pi)
        phase -= 2 * np.pi * (phase > np.pi)
        return phase, hemo


def phase_maps(data, offset_ring=0, offset_wedge=0, do_wedge=True, 
               do_ring=True, do_phase_unwrapping=False, mesh=None, mask=None):
    """ Compute the phase for each functional map
    
    Parameters
    ----------
    data: dictionary with keys 'sin_wedge_pos', 'sin_wedge_neg', 
          'cos_wedge_neg', 'cos_ring_pos', 'sin_ring_neg', 'cos_wedge_pos', 
          'sin_ring_pos', 'cos_ring_neg'
         arrays of shape (n_nodes) showing fMRI activations 
         for different retino conditions
    offset_ring: float, offset value to apply to the ring phase
    offset_wedge: float, offset value to apply to the wedge phase
    mesh: path or mesh instance, underlying mesh model
    mask: array of shape (n_nodes), where to do the analysis 
    do_wedge: bool, should we do the ring phase estimation or not
    do_ring: bool, should we do the ring phase estimation or not
    do_phase_unwrapping: bool, whether or not to correct for 2pi errors
    """
    phase_ring, whase_wedge, hemo = None, None, None
    if do_ring:
        phase_ring_pos = np.arctan2(data['sin_ring_pos'], data['cos_ring_pos'])
        phase_ring_neg = np.arctan2(data['sin_ring_neg'], data['cos_ring_neg'])
        phase_ring, hemo_ring = combine_phase(
            phase_ring_pos, phase_ring_neg, offset_ring, hemo=.6)
        hemo = hemo_ring

    if do_wedge:
        phase_wedge_pos = np.arctan2(data['sin_wedge_pos'], 
                                     data['cos_wedge_pos'])
        phase_wedge_neg = np.arctan2(data['sin_wedge_neg'], 
                                     data['cos_wedge_neg'])
        phase_wedge, hemo_wedge = combine_phase(
            phase_wedge_pos, phase_wedge_neg, offset_wedge, hemo=.6)
        hemo = hemo_wedge
        if do_phase_unwrapping:
            phase_unwrapping(phase_wedge, mesh, mask)

    if do_ring and do_wedge:
        hemo = 0.5 * (hemo_ring + hemo_wedge)

    return phase_wedge, phase_ring, hemo


def phase_unwrapping(phase, mesh, mask, path=None):
    """ This aims at solging 2-pi phase shifts
    """
    for iteri in range(5):
        mg = mep.mesh_to_graph(mesh).subgraph(mask.ravel())
        diff = np.ravel(phase[mg.edges.T[0]] - 
                        phase[mg.edges.T[1]])
        pos_vertex = np.zeros(mg.V)
        pos_vertex[np.unique(mg.edges.T[0, diff > np.pi])] = 1
        pos_vertex[np.unique(mg.edges.T[0, diff < - np.pi])] = - 1
        mg.remove_edges(np.abs(diff) < np.pi / 2)
        u = mg.cc()
        change = False
        for k in np.unique(u):
            if np.sum(u == k) < np.size(u) / 2:
                disc = pos_vertex[u == k].sum() / (
                    1.e-6 + np.abs(pos_vertex[u == k]).sum())
                if np.abs(disc) > .75:
                    change = True
                if (disc > .75):
                    phase[u == k] -= 2 * np.pi   
                if (disc < -.75):
                    phase[u == k] += 2 * np.pi
        if path is not None: 
            ex = - 1 * np.ones(np.size(mask))
            ex[mask.squeeze()] = u
            save_texture(path, ex, 'estimate')
        if change == False:
            break
    phase[phase < - np.pi] = - np.pi + 1.e-6
    phase[phase > np.pi] = np.pi - 1.e-6
    return phase


def visual_field_scalar(mesh, ring, wedge, mask=None):
    """ Compute the visual field scalar at each point of the mesh
    
    Parameters
    ---------
    mesh: string, a path to a mesh file
    ring: string, a path to the ring texture for that mesh/subject
    wedge: string, a path to the wedge texture for that mesh/subject
    mask: string, a path to mask texture for that mesh/subject,
          so that only part of the vertices are conssidered.

    Returns
    -------
    gradient: array of shape (mesh.n_vertices) 
              the visual field sign at the mesh vertices
    """
    normals = mep.compute_normal_vertex(mesh, mask)
    ring_grad = mep.texture_gradient(mesh, ring, mask)
    wedge_grad = mep.texture_gradient(mesh, wedge, mask)
    vfs = np.array([np.linalg.det(np.vstack((g, r, w))) for (g, r, w) in zip (
                normals, ring_grad, wedge_grad)])
    return vfs


def visual_areas_delineation(mesh, ring, wedge, lat=None, lon=None, side=None, 
                             mask=None, vfs=None, areas=None, sigma=3.,
                             default_val=NEGINF):
    """ Try to identify the visual areas using the vfs technique
    
    Parameters
    ---------
    mesh: string, a path to a mesh file
    ring: string, a path to the ring texture for that mesh/subject
    wedge: string, a path to the wedge texture for that mesh/subject
    lat: string, path to the latitude texture
    lon: string, path to the longitude texture
    side: string to be chosen among 'left' or 'right' or None, optional,
          On which hemisphere we are working
    mask: string or None, optional,
        a path to mask texture for that mesh/subject,
        so that only part of the vertices are conssidered.
    vfs: string or None, optional, 
         path where the visual areas will be written
    areas: string or None , optional, 
         path where the visual areas will be written
    sigma: float, optional, size of the smoothing kernel
    default_val: float optional, default value for masked out regions.
    
    Returns
    -------
    v_areas: an array yielding the following labels:
            -1: regions with little retinotopic activity
            0: V1 dorsal
            1: V1 ventral
            2: V2 dorsal
            3: V2 ventral
            4: V3 dorsal
            5: V3 ventral
            6: other retinotopic regions
    """
    # finnesse the mask
    
    # get the visual field sign
    # swedge = op.join(op.dirname(wedge), 's' + op.basename(wedge))
    #_ = mep.smooth_texture(mesh, wedge, output_texture=swedge, sigma=sigma, 
    #                       mask=mask)
    vfd = visual_field_scalar(mesh, ring, wedge, mask)

    # write fs as a texture
    n_nodes = load_texture(ring).size
    if mask == None:
        mask = np.ones(n_nodes).astype(np.bool)
    else: 
        mask = load_texture(mask).astype(np.bool)
    wvfd = 0 * np.ones(n_nodes)
    wvfd[mask] = vfd
    save_texture(vfs, wvfd)
    # get the connected components
    # get their latitude, longitude
    # identify them


def get_ring_minimum(xy, ring, thresh=0):
    """Get the minimum of the ring values, assuming an ellispoidal shape"""
    x, y = xy[ring < thresh].mean(0) 
    print x, y
    """
    from sklearn.covariance.outlier_detection import EllipticEnvelope
    inliers = EllipticEnvelope(contamination=.1).fit(xy).support_ *\
              (ring < thresh)
    xy = xy.copy()[inliers]
    xy2 = np.vstack((np.ones(xy.shape[0]), xy.T[0], xy.T[1], xy.T[0] ** 2, 
                     xy.T[1] ** 2, xy.T[0] * xy.T[1])).T
    A = np.dot(np.linalg.pinv(xy2), ring[inliers])  
    grad = np.array([A[1], A[2]])
    hessian = 2 * np.array([[ A[3], A[5] / 2], [A[5] / 2, A[4]]])
    x, y = - np.dot(np.linalg.inv(hessian), grad)
    print x, y
    """
    return np.array([x, y])


def retino_template(xy, ring, wedge, mesh, mask_, verbose=True, side='left'):
    """"""
    center = get_ring_minimum(xy, ring, 0) # - np.pi / 2)
    ecc = np.sqrt(np.sum((xy - center) ** 2, 1))
    angle = - np.arctan2(xy.T[1] - center[1], xy.T[0] - center[0])
    radius = np.sqrt(np.mean(ecc ** 2)) * 1.2

    def retino_polar(angle, ecc, radius, scale=1.):
        """"""
        a1, a2, a3 = scale * np.pi *.2, scale * np.pi *.3, scale * np.pi * .4
        delta = .2 * 4. / radius ** 2
        corr = angle * (radius ** 2 / 4 - (ecc - radius / 2) ** 2)
        angle_ = angle - delta * corr
        mask = (ecc < radius) * (np.abs(angle_) < a3)
        polar = np.zeros_like(angle_)
        slope = np.pi / (2 * a1)
        b1 = 2 * slope * a1
        polar = slope * angle_
        polar[angle_ > a1] = b1 - slope * angle_[angle_ >  a1]
        polar[angle_ > a2] = slope * angle_[angle_ > a2] - 2 * slope * (a2 - a1)
        
        polar[angle_ < - a1] = - b1 - slope * angle_[angle_ < - a1]
        polar[angle_ < - a2] = slope * (angle_[angle_ < - a2] + 2 * (a2 - a1))
        polar[mask == 0] = 0
        maps = {'V1v': mask * (angle_ < 0) * (angle_ > - a1),
                'V1d': mask * (angle_ > 0) * (angle_ < a1),
                'V2v': mask * (angle_ < -a1) * (angle_ > - a2),
                'V2d': mask * (angle_ > a1) * (angle_ < a2),
                'V3v': mask * (angle_ < - a2) * (angle_ > - a3),
                'V3d': mask * (angle_ > a2) * (angle_ < a3),
                }
        return mask, polar, maps

    vmin, vmax = 0., np.pi
    if side == 'left':
        vmin, vmax = -np.pi, 0.

    best_score = - np.inf
    # 2-dimensional search of the best rotation/scaling
    for scale in [1.]:
        for theta in np.linspace(0, 2 * np.pi, 100):
            angle_ = angle + theta
            angle_[angle_ > np.pi] -= 2 * np.pi
            mask, polar, maps = retino_polar(angle_, ecc, radius, scale)
            if mask.sum() == 0:
                continue
            if side == 'left':
                polar -= np.pi / 2
            else:
                polar += np.pi / 2
            weight = mask * (polar > vmin) * (polar < vmax)
            score =  1. - np.sum(
                weight[mask] * (wedge[mask] - polar[mask]) ** 2)/\
                np.sum(weight[mask] * (wedge[mask]) ** 2)
            score = np.corrcoef(wedge[mask], polar[mask])[0, 1]
            if score > best_score:
                best_score = score
                best_polar = polar
                best_theta = theta
                visual_maps = maps
                best_mask = mask
                best_corr = np.corrcoef(wedge[mask], polar[mask])[0, 1]
    
    """
    # finer search
    for scale in [0.6, 0.8, 1.]:
        for theta in np.linspace(best_theta - .5, best_theta + .5, 21):
            angle_ = angle + theta
            angle_[angle_ > np.pi] -= 2 * np.pi
            angle_[angle_ > np.pi] -= 2 * np.pi
            angle_[angle_ < - np.pi] += 2 * np.pi
            mask, polar, maps = retino_polar(angle_, ecc, radius, scale)
            if mask.sum() == 0:
                continue
            if side == 'left':
                polar -= np.pi / 2
            else:
                polar += np.pi / 2
            weight = mask * (polar > vmin) * (polar < vmax)
            score =  1. - np.sum(
                weight[mask] * (wedge[mask] - polar[mask]) ** 2) / np.sum(
                weight[mask] * (wedge[mask]) ** 2)
            if score > best_score:
                best_score = score
                best_polar = polar
                best_theta = theta
                visual_maps = maps
                best_mask = mask
                best_corr = np.corrcoef(wedge[mask], polar[mask])[0, 1]
    """

    if verbose:
        print best_theta, best_corr    
        import matplotlib.pyplot as plt
        plt.figure(figsize = (10, 5))
        coord, tri = mep.mesh_arrays(mesh)
        coord, tri = coord[mask_], tri[mask_[tri].all(1)]
        tri = np.hstack((0, np.cumsum(mask_)[:-1]))[tri]
        plt.subplot(1, 3, 1)
        plt.tripcolor(xy.T[0], xy.T[1], tri, best_polar * best_mask, 
                      vmin=vmin, vmax=vmax)
        plt.plot(center[0], center[1], '+k', linewidth=4)
        plt.subplot(1, 3, 2)
        plt.tripcolor(xy.T[0], xy.T[1], tri, wedge, vmin=vmin, vmax=vmax)
        plt.plot(center[0], center[1], '+k', linewidth=4)
        plt.subplot(1, 3, 3)
        plt.tripcolor(xy.T[0], xy.T[1], tri, ring)
        plt.plot(center[0], center[1], '+k', linewidth=4)
    
    return mask, polar, visual_maps


##################################################################
# Ancillary functions
##################################################################


def angular_maps(side, paths, all_reg, threshold=3.1, size_threshold=10,
                 offset_wedge=0, offset_ring=0, smooth=0, 
                 do_wedge=True, do_ring=True, do_phase_unwrapping=False):
    """
    Parameters
    ----------
    side: string, either 'left', 'right' or False
    paths: dictionary, set of paths generated during the GLM analysis
    all_reg: list of strings,
            identifiers of the contrast files used in angular mapping
    threshold: float, optional
               threshold defining the brain regions
               where the analysis is performed
    size_threshold: int, optional
                    threshold on connected components size
                    to remove isolated voxels
    offset_wedge: float, optional
                  offset to be applied to wedge angle
    offset_ring float, optional
                  offset to be applied to ring angle
    """
    #-------------------------------------------------------------------------
    # First get the region with significant activity at the stimuli frequency
    #-------------------------------------------------------------------------

    if side == False:
        stat_map = op.join(paths["contrasts"],
                                'effects_of_interest_z_map.nii')
 
        ## create an occipital data_mask
        mask = load(stat_map).get_data() > threshold

        # remove isolated voxels
        mask = cc_array_mask(mask, size_threshold)

        # load and mask the data
        data = {}
        for r in all_reg:
            contrast_file = op.join(paths["contrasts"],
                                             '%s_con.nii' % r)
            data[r] = load(contrast_file).get_data()[mask]
        do_phase_unwrapping = False
        mesh = None
    else:
        const_map = op.join(paths["contrasts"],
                                 '%s_constant_con.gii' % side)
        if op.exists(const_map):
            const = load_texture(const_map)
        else:
            const = 0
        mesh = paths['%s_mesh' % side]
        stat_map = op.join(paths["contrasts"],
                                '%s_effects_of_interest_z_map.gii' % side)
        smooth_stat_map = op.join(paths["contrasts"],
                                       '%s_effects_of_interest_z_map_smooth.gii'
                                       % side)
        if smooth > 0:
            stat_map = mep.smooth_texture(
                mesh, stat_map, smooth_stat_map, smooth)
        else:
            stat_map = load_texture(stat_map)
        mask = (np.ravel(stat_map) > threshold) * (np.isnan(const) == 0)
        mask = cc_mesh_mask(mesh, mask, size_threshold).ravel()

        data = {}
        for r in all_reg:
            contrast_file = op.join(
                paths["contrasts"], '%s_%s_con.gii' % (side, r))
            if smooth == 0:
                data[r] = load_texture(contrast_file).ravel()[mask.ravel()]
            else:
                data[r] = mep.smooth_texture(
                    mesh, contrast_file, None, smooth)[mask.ravel()]
                
    # Then compute the activation phase in these regions
    phase_wedge, phase_ring, hemo = phase_maps(
        data, offset_ring, offset_wedge, do_wedge, do_ring, 
        do_phase_unwrapping, mesh=mesh, mask=mask)

    # delineate the visual areas
    if side != False:
        planar_coord = mep.isomap_patch(mesh, mask)
        visual_areas = delineate_areas(phase_ring, phase_wedge,  hemo, mesh, 
                                       mask, planar_coord, side)
  
    # write the results
    data_, id_ = [hemo, mask[mask > 0]], ['hemo', 'mask']
    if do_ring:
        data_.append(phase_ring)
        id_.append('phase_ring')
    if do_wedge:
        data_.append(phase_wedge)
        id_.append('phase_wedge')
    if side != False:
        data_.append(visual_areas)
        id_.append('visual_areas')
        
    if side == False:
        for (x, name) in zip(data_, id_):
            wdata = np.zeros(load(stat_map).shape)
            wdata[mask > 0] = x
            wim = Nifti1Image(wdata, load(stat_map).get_affine())
            save(wim, op.join(paths["contrasts"], '%s.nii' % name))
    else:
        for (x, name) in zip(data_, id_):
            ex = 0 * np.ones(np.size(mask)).astype(np.float32)
            ex[mask.squeeze()] = x
            save_texture(
                op.join(paths["contrasts"], side + '_' + name + '.gii'),
                ex, 'estimate')


def delineate_areas(phase_ring, phase_wedge, hemo, mesh, mask, planar_coord, 
                    side, write_dir=None):
    """Delineate visual areas

    Parameters
    ---------
    phase_ring: array of shape(n_nodes): phase values for the ring condition
    phase_wedge: array of shape(n_nodes): phase values for the wedge condition
    hemo: array of shape(n_nodes): phase values for the hemodynamic delay
    mask: array of shape (n_nodes), where to do the analysis 
    mesh: path or mesh instance, underlying mesh model
    planar_coord:  array of shape(n_nodes, 2), coordinate system on the surface 
    write_dir: string, where to write the result
    side: string 'left' or 'right'
    """
    _, polar, visual_maps = retino_template(
        planar_coord, phase_ring, phase_wedge, mesh, mask, side=side)
    visual_areas = np.zeros_like(mask).astype(np.float32)
    for idx, key in enumerate(['V1v', 'V1d', 'V2v', 'V2d', 'V3v', 'V3d']):
        visual_areas[mask] += (idx + 1) * visual_maps[key]

    if write_dir is not None:
        # write the results
        phase_wedge = phase_wedge.astype(np.float32)
        phase_ring = phase_ring.astype(np.float32)
        data_, id_ = ([mask[mask > 0], hemo, phase_ring, phase_wedge, 
                       visual_areas[mask]],
                      ['mask', 'hemo', 'phase_ring', 'phase_wedge', 'areas'])
        for (x, name) in zip(data_, id_):
            ex = 0 * np.ones(np.size(mask)).astype(np.float32)
            ex[mask.squeeze()] = x
            save_texture(
                op.join(write_dir, side + '_' + name + '.gii'), ex, 'estimate')

    return visual_areas
