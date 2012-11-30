# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
This script eprforms the same hings as script_glm_look_loc.py
But it uses a block model for the glm

Author: Bertrand Thirion, 2010-2012
"""

import os
import numpy as np
import pylab
import glob

from nibabel import load, save, Nifti1Image

from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
from nipy.modalities.fmri.glm import GeneralLinearModel, data_scaling
from nipy.labs import compute_mask_files

from retino.angular_analysis import (
    load_texture, save_texture, cc_array_mask, cc_mesh_mask, phase_unwrapping)
from config_look_loc import make_paths

# -----------------------------------------------------------
# --------- Set the paths -----------------------------------
#-----------------------------------------------------------

paths = make_paths()
subjects = paths.keys()
result_dir = 'block'

# choose volume-based or surface-based analysis
sides = ['left', 'right'] # [False] #
# False: volume-based analysis
# left: left hemisphere
# right: right hemisphere

# ---------------------------------------------------------
# -------- General data-related Information ---------------
# ---------------------------------------------------------
tr = 2.4
nb_frames = 166
frametimes = np.arange(nb_frames) * tr
n_models = 30
duration, length = 10., 387.

def make_paradigm(delay, duration=10., length=387., odd=False):
    """ Writes the paradigm file for a certain delay"""
    period = 30 # in seconds
    if odd:
        delay = period - delay
    onsets = np.hstack((delay + period * np.arange(6), length - delay
                        - period * np.arange(6) - duration))
    condition_id = ['checkerboard']
    cids = condition_id * len(onsets)
    durations = duration * np.ones(len(onsets))
    return BlockParadigm(cids, onsets, durations)

contrast_val = np.eye(8)[0]

# ---------------------------------------------------------
# ------ First level analysis parameters ---------------------
# ---------------------------------------------------------

# Masking parameters
inf_threshold, sup_threshold = 0.7, .9

# Design matrix
drift_model = "Cosine"
hf_cut = 128

# Treat sequentially all subjects & acquisitions
for subject in subjects[:1]:
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
        
        for (sess, (session, fmri_data)) in enumerate(zip(
            sessions, fmri_series)):
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
                
            best_z, best_d = - np.ones(Y.shape[1]), - np.ones(Y.shape[1])
            
            # create design matrices
            for i_delay, delay in enumerate(
                np.linspace(0, 30, n_models + 1)[:-1]):
                paradigm = make_paradigm(
                    delay, duration=duration, length=length, 
                    odd=(subject in ['sujet10', 'sujet12']))
                design_matrix = make_dmtx(
                    frametimes, paradigm, drift_model=drift_model, hfcut=hf_cut)
                
                # plot the design matrix
                ax = design_matrix.show()
                ax.set_position([.05, .25, .9, .65])
                ax.set_title('Design matrix')

                # fit the glm
                print 'Fitting a GLM (this takes time)...'
                result = GeneralLinearModel(design_matrix.matrix)
                result.fit(Y, model='ar1', steps=100)
                z_values = result.contrast([contrast_val]).z_score()
                best_d[z_values > best_z] = i_delay 
                best_z = np.maximum(z_values, best_z)

            write_array = mask_array.astype(np.float)
            write_array[mask_array] = best_z
            best_d = best_d * 2 * np.pi / n_models
            best_d[best_d > np.pi] -= 2 * np.pi
            if side == False:
                contrast_path = os.path.join(write_dir, '%s_z_map.nii' %
                                             session)
                save(Nifti1Image(write_array, affine), contrast_path)
                write_array[mask_array] = best_d
                contrast_path = os.path.join(write_dir, '%s_d_map.nii' %
                                             session)
                save(Nifti1Image(write_array, affine), contrast_path)
            else:
                contrast_path = os.path.join(write_dir, '%s_%s_z_map.gii' %
                                             (side, session))
                save_texture(contrast_path, write_array)
                contrast_path = os.path.join(write_dir, '%s_%s_d_map.gii' %
                                             (side, session))
                write_array[mask_array] = best_d.astype(np.float)
                save_texture(contrast_path, write_array)
            print (session, duration, length, np.sum(best_z > 3.),
                   np.sum(best_z[best_z > 3.]), np.mean(best_z[best_z > 3.]))


#----------------------------------------------------------------
# Retinotopy specific analysis: phase maps
#----------------------------------------------------------------
import os.path as op


for subject in subjects:
    print ('Computing phase maps in subject %s' % subject)
    # offset_wedge and offset_ring are related to the actual stimulus
    offset_wedge = - np.pi / 2
    # means that the wedge starts at the right horizontal position
    offset_ring = np.pi
    # means that the ring starts at the middle position
    threshold = 3.0
    size_threshold = 50
    do_phase_unwrapping = True
    wedge_path = [p for p in paths[subject]['sessions'] if 'wedge' in p][0]
    write_dir = os.path.join(paths[subject]['base'], 
                            paths[subject]['acquisition'], result_dir)
    if side == False:
        ## create an occipital data_mask
        stat_map = op.join(write_dir, '%s_%s_z_map.nii'  % wedge_path)
        mask_array = load(stat_map).get_data() > threshold
        mask_array = cc_array_mask(mask_array, size_threshold)
        #
        # load and mask the data
        contrast_file = op.join(write_dir, '%s_d_map.nii' % wedge_path)
        data = load(contrast_file).get_data()[mask_array] + (
            offset_wedge + np.pi / 2 + np.pi / 12)
        data[data > np.pi] -= 2 * np.pi
        data[data < -np.pi] += 2  * np.pi
        write_array = mask_array.astype(np.float)
        write_array[mask_array] = data
        contrast_path = os.path.join(write_dir, 'phase_wedge.nii')
        save(Nifti1Image(write_array, affine), contrast_path)
        sides = []
    #
    for side in sides:
        mesh = paths[subject]['%s_mesh' % side]
        stat_map = op.join(write_dir, '%s_%s_z_map.gii' % (side, wedge_path))
        contrast_map = op.join(write_dir, '%s_%s_d_map.gii' % (side,
                                                               wedge_path))
        stat_map = load_texture(stat_map)
        mask_array = cc_mesh_mask(mesh, np.ravel(stat_map) > threshold,
                            size_threshold)
        data = load_texture(contrast_map).ravel()[mask_array] + (
            offset_wedge +  np.pi / 2 + np.pi / 12)
        data[data > np.pi] -= 2 * np.pi
        data[data < -np.pi] += 2  * np.pi
        if do_phase_unwrapping:
            phase_unwrapping(data, mesh, mask_array)
        write_array = mask_array.astype(np.float)
        write_array[mask_array] = data
        contrast_path = os.path.join(write_dir, '%s_phase_wedge.gii' % side)
        save_texture(contrast_path, write_array)
        write_array = mask_array.astype(np.float)
        contrast_path = os.path.join(write_dir, '%s_mask.gii' % side)
        save_texture(contrast_path, write_array)

pylab.show()
