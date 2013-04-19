# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Script that performs the GLM analysis on the cortical surface
In order to obtain activation  maps and then retinotopic maps that yield
activation phase-based information.
It relies on the nipy library.

Author: Bertrand Thirion, 2010-2012
"""

import os
import numpy as np
import pylab
import glob

from nibabel import load, save, Nifti1Image
from retino.angular_analysis import load_texture, save_texture, angular_maps

from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.glm import GeneralLinearModel, data_scaling
from nipy.labs import compute_mask_files
from config_look_loc import make_paths

# -----------------------------------------------------------
# --------- Set the paths -----------------------------------
#-----------------------------------------------------------

paths = make_paths()
subjects = paths.keys()[:1]
result_dir = 'analysis'

# choose volume-based or surface-based analysis
sides = ['left', 'right'] #[False] #
# False: volume-based analysis
# left: left hemisphere
# right: right hemisphere

# ---------------------------------------------------------
# -------- General data-related Information ---------------
# ---------------------------------------------------------
tr = 2.4
nb_frames = 166
frametimes = np.arange(nb_frames) * tr

period = 30 # in seconds
r1 = np.sin(2 * np.pi * 1. / period * frametimes) * (frametimes < 180)
r2 = np.cos(2 * np.pi * 1. / period * frametimes) * (frametimes < 180)
r3 = np.sin(2 * np.pi * 1. / period * frametimes) * (frametimes > 210)
r4 = np.cos(2 * np.pi * 1. / period * frametimes) * (frametimes > 210)

reg_matrix = np.vstack((r1, r2, r3, r4)).T

all_reg = ['sin_ring_pos', 'cos_ring_pos', 'sin_ring_neg', 'cos_ring_neg',
           'sin_wedge_neg', 'cos_wedge_neg', 'sin_wedge_pos', 'cos_wedge_pos']

# ---------------------------------------------------------
# ------ First level analysis parameters ---------------------
# ---------------------------------------------------------

# Masking parameters
inf_threshold, sup_threshold = 0.7, .9

# Design matrix
drift_model = "Cosine"
hf_cut = 128


def make_contrasts(sessions, n_reg=11, odd=False):
    """ Build the contrasts for the loolkloc experiment

    Parameters
    ==========
    sessions: list of strings, ids of the sessions
    n_reg: int, number of regressors in the deisgn matrix
    odd: boolean, 1 if positive first

    It is assumed that the sessions ids correspond to the
    (ring session, wedge session)
    """
    ring, wedge = sessions
    contrasts = {
        'ring': {ring: np.eye(n_reg)[: 4], wedge: np.zeros((n_reg, 4))},
        'wedge': {wedge: np.eye(n_reg)[: 4], ring: np.zeros((n_reg, 4))},
        'effects_of_interest': {
            ring: np.eye(n_reg)[:4],
            wedge: np.eye(n_reg)[:4]},
        'sin_wedge_pos': {ring: np.zeros(n_reg), wedge: np.eye(n_reg)[0]},
        'sin_wedge_neg': {ring: np.zeros(n_reg), wedge: np.eye(n_reg)[2]},
        'cos_wedge_pos': {ring: np.zeros(n_reg), wedge: np.eye(n_reg)[1]},
        'cos_wedge_neg': {ring: np.zeros(n_reg), wedge: np.eye(n_reg)[3]},
        'sin_ring_pos': {wedge: np.zeros(n_reg), ring: np.eye(n_reg)[0]},
        'sin_ring_neg': {wedge: np.zeros(n_reg), ring: np.eye(n_reg)[2]},
        'cos_ring_pos': {wedge: np.zeros(n_reg), ring: np.eye(n_reg)[1]},
        'cos_ring_neg': {wedge: np.zeros(n_reg), ring: np.eye(n_reg)[3]},
        'constant': {ring: np.eye(n_reg)[-1], wedge: np.eye(n_reg)[-1]},
        }
    if odd:
        contrasts['sin_wedge_pos'][wedge] = np.eye(n_reg)[2]
        contrasts['sin_wedge_neg'][wedge] = np.eye(n_reg)[0]
        contrasts['cos_wedge_pos'][wedge] = np.eye(n_reg)[3]
        contrasts['cos_wedge_neg'][wedge] = np.eye(n_reg)[1]
    return contrasts


# Treat sequentially all subjects & acquisitions
for subject in subjects:
    #  generate all paths
    fmri_dir = os.path.join(paths[subject]['base'],
                            paths[subject]['acquisition'])
    write_dir = os.path.join(fmri_dir, result_dir)
    epi_mask = os.path.join(write_dir, 'mask.nii')
    sessions = paths[subject]['sessions']

    if not os.path.exists(write_dir):
        os.mkdir(write_dir)

    # main loop on sides
    for side in sides:
        wild_card = paths[subject]['wild_card']
        if side == 'left':
            wild_card = '*lh_smooth5.gii' # '*lh.gii' #
        elif side == 'right':
            wild_card = '*rh_smooth5.gii' # '*rh.gii' #

        print "Subject : %s, side : %s" % (subject, side)
        # get the images
        fmri_series = [glob.glob(
                os.path.join(fmri_dir, session, wild_card))
                       for session in sessions]

        # warning: glob returns the files in random order
        for f in fmri_series:
            f.sort()

        if side == False:
            # in the volume, compute a mask o fthe brain
            mask_array = compute_mask_files(fmri_series[0], epi_mask, True,
                                            inf_threshold, sup_threshold)[0]

        # Specify the contrasts
        contrasts = make_contrasts(sessions, odd=subject in
                                   ['sujet10', 'sujet12'])

        contrast_obj = {}
        for (sess, (session, fmri_data)) in enumerate(zip(
            sessions, fmri_series)):
            # create design matrices
            reg = all_reg[4 * sess: 4 * sess + 4] # fixme
            design_matrix = make_dmtx(
                frametimes, add_regs=reg_matrix, add_reg_names=reg,
                drift_model=drift_model, hfcut=hf_cut)

            # plot the design matrix
            ax = design_matrix.show()
            ax.set_position([.05, .25, .9, .65])
            ax.set_title('Design matrix')
            pylab.savefig(os.path.join(write_dir, 'design_matrix_%s.png') %\
                              session)
            # get the data
            if side == False:
                Y, _ = data_scaling(np.array([load(f).get_data()[mask_array]
                                              for f in fmri_data]))
                affine = load(fmri_data[0]).get_affine()
            else:
                Y, _ = data_scaling(np.array(
                        [load_texture(f) for f in fmri_data]))
                mask_array = np.var(Y, 0) > 0
                Y = Y[:, mask_array]

            # fit the glm
            print 'Fitting a GLM (this takes time)...'
            result = GeneralLinearModel(design_matrix.matrix)
            result.fit(Y, model='ar1', steps=100)

            for contrast_id, contrast_val in contrasts.iteritems():
                if (contrast_val[session] == 0).all():
                    continue

                if contrast_id in contrast_obj.keys():
                    contrast_obj[contrast_id] = contrast_obj[contrast_id] +\
                        result.contrast(contrast_val[session])
                else:
                    contrast_obj[contrast_id] =\
                        result.contrast(contrast_val[session])

        print 'Computing contrasts...'
        # In fact, this only writes the output images/textures
        for index, contrast_id in enumerate(contrasts):
            if contrast_id == 'effects_of_interest':
                contrast_obj[contrast_id].stat_ = 0.5 * (
                    contrast_obj['wedge'].stat_ + contrast_obj['ring'].stat_)
                contrast_obj[contrast_id].p_value_ = None # recompute it
            z_values = contrast_obj[contrast_id].z_score()
            effect_ = contrast_obj[contrast_id].effect
            print '  Contrast % 2i out of %i: %s' % (
                index + 1, len(contrasts), contrast_id)

            if side == False:
                contrast_path = os.path.join(write_dir, '%s_z_map.nii' %
                                             contrast_id)
                write_array = mask_array.astype(np.float)
                write_array[mask_array] = z_values
                contrast_image = Nifti1Image(write_array, affine)
                save(contrast_image, contrast_path)

                if effect_.shape[0] == 1:
                    write_array[mask_array] = np.ravel(effect_)
                    contrast_image = Nifti1Image(write_array, affine)
                    save(contrast_image, os.path.join(write_dir, '%s_con.nii' %
                                                      contrast_id))
            else:
                contrast_path = os.path.join(write_dir, '%s_%s_z_map.gii' %
                                             (side, contrast_id))
                write_array = mask_array.astype(np.float)
                write_array[mask_array] = z_values
                save_texture(contrast_path, write_array)
                contrast_path = os.path.join(write_dir, '%s_%s_con.gii' %
                                             (side, contrast_id))
                if effect_.shape[0] == 1:
                    write_array[mask_array] = np.ravel(effect_)
                    save_texture(contrast_path, write_array)

#----------------------------------------------------------------
# Retinotopy specific analysis: phase maps
#----------------------------------------------------------------

for subject in subjects:
    if sides == False:
        continue
    print ('Computing phase maps in subject %s' % subject)
    contrast_path = os.path.join(
            paths[subject]['base'], paths[subject]['acquisition'], result_dir)
    #
    # offset_wedge and offset_ring are related to the actual stimulus
    offset_wedge = - np.pi / 2
    # means that the wedge starts at the right horizontal position
    offset_ring = np.pi
    # means that the ring starts at the middle position
    threshold = 3.0

    for side in sides:
        if side is not False:
            mesh_path = make_paths()[subject]['%s_mesh' % side]
        # threshold of the main effects map, in z-value
        angular_maps(
            side, contrast_path, mesh_path, 
            all_reg=all_reg, 
            threshold=threshold,
            size_threshold=50, 
            offset_wedge=offset_wedge,
            offset_ring=offset_ring, 
            smooth=0., 
            do_phase_unwrapping=True)
