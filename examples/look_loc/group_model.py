"""
Make a group analysis out of the individual data

Note
====
first type export SUBJECTS_DIR=/volatile/thirion/data/lookloc/

Todo
====

Author:Bertrand Thirion, 2012
"""

import commands
import numpy as np
import os.path as op
from nibabel.freesurfer import mghformat
from nibabel import Nifti1Image
import enthought.mayavi.mlab as mlab
from parietal.glm_files_layout.cortical_glm import load_texture
from config_look_loc import make_paths
from angular_analysis import main_cc, phase_maps, delineate_areas
from visu_mlab import plot_retino_image
try:
    from parietal.surface_operations.mesh_processing import isomap_patch
except:
    from parietal_copy.mesh_processing import isomap_patch


def resample_to_average(tex_path, subject, side, verbose=False):
    """Resample a given texture to fsaverage

    Parameters
    ==========
    tex_path: string, path of the input texture
    subject: string, subject id in the freesurfer database
    side: string, on of ['left', 'right']
    verbose: boolean, the verbosity mode
    
    Returns
    =======
    resmapled: string, path of the result
    """
    resampled = op.join(op.dirname(tex_path), 
                        op.basename(tex_path)[:-4] + '_resampled.gii')

    # convert the input to .mgz format
    mgz = tex_path[:-4] + '.mgz'
    tex = load_texture(tex_path)[np.newaxis, np.newaxis].T
    mghformat.save(Nifti1Image(tex, np.eye(4)), mgz)

    # run the resmapling using freesurfer tools
    fs_comment = commands.getoutput(
        'mri_surf2surf --srcsubject %s --srcsurfval %s --trgsubject ico --trgsurfval %s --hemi %sh --trgicoorder 7' % (
            subject, mgz, resampled, side[0]))
    if verbose:
        print fs_comment
    return resampled


def resample_to_subject(tex_path, subject, side, output_path, verbose=False):
    """Resample a given texture from fsaverage to subject space
    
    Parameters
    ==========
    tex_path: string, path of the input texture
    subject: string, subject id in the freesurfer database
    side: string, on of ['left', 'right']
    outoput_path: string
    verbose: boolean, the verbosity mode
    """
    # convert the input to .mgz format
    mgz = tex_path[:-4] + '.mgz'
    tex = load_texture(tex_path)[np.newaxis, np.newaxis].T
    mghformat.save(Nifti1Image(tex, np.eye(4)), mgz)
    fs_comment = commands.getoutput(
        'mri_surf2surf --trgsubject %s --trgsurfval %s --srcsubject ico --srcsurfval %s --hemi %sh --srcicoorder 7' % (
                subject, output_path, mgz, side[0]))
    if verbose:
        print fs_comment
    return output_path


def register_and_average_data(
    data=None, positions=None, save=False, path='/tmp/data.npz', 
    cor_path='/tmp/corrected_data.npz', n_samples=1000):
    """Performs some registration of the data using a daemon algorithm"""
    from parietal.registration.daemons import demons_algo, learn_scale
    from sklearn.metrics.pairwise import rbf_kernel

    if save:
        np.savez(path, data=data, positions=positions)
    if data == None:
        data = np.load(path)['data']
    if positions == None:
        positions = np.load(path)['positions']

    def equalize_histo(source, template):
        """Equalizes the histogram of source to that of template"""
        stemplate = np.sort(template, 0)
        index = np.argsort(source, 0)
        source_ = np.zeros_like(source)
        for i in range(source.shape[1]):
            source_[index[:, i], i] = stemplate[:, i]
        return source_

    keys = data.keys()
    # compute the mean data
    mean_data = dict([(r, np.mean(data[r], 0)) for r in keys])
    
    avg_data = np.array([mean_data[r] for r in keys]).T
    
    # register each subject to the template
    p_range = np.array([.5, .75, 1., 1.25])
    # p_range = np.array([1., 1.25, 1.5, 2.])
    for subject in range(data[r].shape[0]):
        aux = np.random.permutation(data[r].shape[1])[:n_samples]
        scale = learn_scale(positions[aux])
        smat = rbf_kernel(positions, gamma=1. / scale ** 2)
        smat /= smat.sum(0)
        subject_data = np.array([np.dot(smat.T, data[r][subject]) 
                                 for r in keys]).T
        subject_data = equalize_histo(subject_data, avg_data)
        """#
        disp, conv, param_ = demons_algo(
            positions[aux], subject_data[aux], avg_data[aux], scale=scale, 
            param_range=p_range, n_iter=20, delta=1.e-4)
        print subject, param_
        subject_data = conv.predict(positions)
       #
        for i, r in enumerate(keys):
            data[r][subject] = subject_data[:, i]    

    """# recompute the mean_data
    mean_data = dict([(r, np.mean(data[r], 0)) for r in keys])
    return mean_data


#########################################################################
# Main script
#########################################################################
paths = make_paths()
subjects = paths.keys()
all_reg = ['sin_ring_pos', 'cos_ring_pos', 'sin_ring_neg', 'cos_ring_neg',
           'sin_wedge_neg', 'cos_wedge_neg', 'sin_wedge_pos', 'cos_wedge_pos']
write_dir = '/volatile/thirion/data/lookloc/fsaverage'
threshold = 10.0
size_threshold = 50
offset_ring, offset_wedge = np.pi, - np.pi / 2
do_phase_unwrapping = True

for side in ['left', 'right']:
    # First compute a mask of the active regions, by thresholding the
    # so-called "average stat" map
    avg_stat = None
    mesh = op.join(op.dirname(paths[subjects[0]]['base']), 'fsaverage',
                   '%sh.inflated.gii' % side[0])
    for subject in subjects:
        tex_dir = op.join(paths[subject]['base'], 
                          paths[subject]['acquisition'], 'analysis')
        stat_map = op.join(tex_dir,
                           '%s_effects_of_interest_z_map.gii' % side)
        r_stat_map = resample_to_average(stat_map, subject, side)
        if avg_stat == None:
            avg_stat = load_texture(r_stat_map)
        else:
            avg_stat += load_texture(r_stat_map)
    mask = avg_stat / np.sqrt(len(subjects)) > threshold
    mask = main_cc(mesh, mask)

    # Then, given the mask, extract the beta data from each subject
    data = dict([ (r, []) for r in all_reg])
    for subject in subjects:
        tex_dir = op.join(paths[subject]['base'], 
                          paths[subject]['acquisition'], 'analysis')
        for r in all_reg:
            contrast_file = op.join(tex_dir, '%s_%s_con.gii' % (side, r))
            r_cfile = resample_to_average(contrast_file, subject, side)
            datar = load_texture(r_cfile)[mask]
            data[r].append(datar)
    
    for r in data.keys():
        data[r] = np.array(data[r])
    
    # Compute the average value of the data
    planar_coord = isomap_patch(mesh, mask)
    data = register_and_average_data(
        data, planar_coord, save=True, path='/tmp/data_%s.npz' % side, 
        cor_path='/tmp/corrected_data_%s.npz' % side)
    phase_wedge, phase_ring, hemo = phase_maps(
        data, offset_ring, offset_wedge, do_wedge=True, do_ring=True, 
        do_phase_unwrapping=do_phase_unwrapping, mesh=mesh, mask=mask)
    visual_areas = delineate_areas(phase_ring, phase_wedge,  hemo, mesh, mask, 
                                   planar_coord, side, write_dir)
    
    #------------------------------------------------------------------------
    # plot the result
    #------------------------------------------------------------------------
    wedge_tex = op.join(write_dir, side + '_phase_wedge.gii')
    curv_tex = op.join(op.dirname(paths[subject]['base']), 'fsaverage',
                   '%sh.curv.gii' % side[0])

    f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
    mlab.clf()
    if side == 'left':
        vmin, vmax = -3.14, .5
    else:
        vmin, vmax = -.5, np.pi
    plot_retino_image(mesh, name="%s_hemisphere" %side, tf=None, 
                      curv = load_texture(curv_tex), mask=mask,
                      tex=load_texture(wedge_tex), vmin=vmin, vmax=vmax)
    mlab.view(280, 120)
    f = mlab.figure(bgcolor=(.05, 0, .1), size=(400, 400))
    mlab.clf()
    plot_retino_image(mesh, name="LeftHemisphere", tf=None, mask=mask,
                      tex=visual_areas, curv=load_texture(curv_tex), vmin=0, 
                      vmax=6)
    mlab.view(280, 120)

    #------------------------------------------------------------------------
    # Resample to individual space
    #------------------------------------------------------------------------
    visual_area_tex = op.join(write_dir, side + '_areas.gii')
    for subject in subjects:
        output_path = op.join(
            paths[subject]['base'], paths[subject]['acquisition'], 'analysis',
            '%s_areas.gii' % side)
        resample_to_subject(visual_area_tex, subject, side, output_path, 
                            verbose=True)
    
