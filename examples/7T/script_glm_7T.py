# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Script that performs the GLM analysis on the cortical surface
In order to obtain activation  maps and then retinotopic maps that yield
activation phase-based information.
It relies quite heavily on nipy library.

fixme
-----
io problems for surface-based analysis

Author: Bertrand Thirion, 2010
"""

import os
import numpy as np
import glob
import pylab
import scipy.stats as st

from nipy.labs import compute_mask_files
from nibabel import load, save, Nifti1Image
try:
    from parietal.glm_files_layout.cortical_glm import (load_texture,
                                                        save_texture)
except:
    from parietal_copy.cortical_glm import load_texture, save_texture
from nipy.modalities.fmri.glm import GeneralLinearModel, data_scaling
from retino import angular_analysis
from config_retino_7T import make_paths
from nipy.modalities.fmri.design_matrix import make_dmtx
from nibabel.gifti import read

# -----------------------------------------------------------
# --------- Set the paths -----------------------------------
#-----------------------------------------------------------

paths = make_paths()
subjects = paths.keys()
result_dir = 'analysis'

# choose volume-based or surface-based analysis
sides =  ['left', 'right']# [False] #
# False: volume-based analysis
# left: left hemisphere
# right: right hemisphere

# generic image name
model_id = 'no_warp_moco'
# wild_card = 'auf*.img'

# ---------------------------------------------------------
# -------- General data-related Information ---------------
# ---------------------------------------------------------

tr = 2.4
nb_frames = 132
frametimes = np.arange(nb_frames) * tr

period = 8.
r1 = np.sin(2 * np.pi * 1. / 38.4 * frametimes)
r2 = np.cos(2 * np.pi * 1. / 38.4 * frametimes)
reg_matrix = np.vstack((r1, r2)).T

# ---------------------------------------------------------
# ------ First level analysis parameters ---------------------
# ---------------------------------------------------------

#---------- Masking parameters
inf_threshold = 0.6
sup_threshold = 0.9

#---------- Design Matrix
drift_model = "Cosine"
hfcut = 128

def make_contrasts(sessions, n_reg=7):
    """ Build the contrasts for the loolkloc experiment

    Parameters
    ==========
    sessions: list of strings, ids of the sessions
    n_reg: int, number of regressors in the deisgn matrix
    It is assumed that the sessions ids correspond to the
    (positive wedge, negative wedge, positive ring, negative ring)
    """
    ring_pos, ring_neg, wedge_pos, wedge_neg = sessions
    con_ids = ['sin_ring_pos', 'cos_ring_pos', 'sin_ring_neg',  'cos_ring_neg',
               'sin_wedge_pos', 'cos_wedge_pos', 'sin_wedge_neg',
               'cos_wedge_neg']

    #Create dummy contrasts
    contrast = dict([(session, np.zeros(n_reg)) for session in sessions])
    contrasts = dict([(con_id, contrast.copy()) for con_id in con_ids])
    contrasts['sin_wedge_pos'][wedge_pos] = np.eye(n_reg)[0]
    contrasts['cos_wedge_pos'][wedge_pos] = np.eye(n_reg)[1]
    contrasts['sin_wedge_neg'][wedge_neg] = np.eye(n_reg)[0]
    contrasts['cos_wedge_neg'][wedge_neg] = np.eye(n_reg)[1]
    contrasts['sin_ring_pos'][ring_pos] = np.eye(n_reg)[0]
    contrasts['cos_ring_pos'][ring_pos] = np.eye(n_reg)[1]
    contrasts['sin_ring_neg'][ring_neg] = np.eye(n_reg)[0]
    contrasts['cos_ring_neg'][ring_neg] = np.eye(n_reg)[1]
    return contrasts, con_ids


# Treat sequentially all subjects & acquisitions
for subject in subjects:
    for side in sides:
        print("Subject: %s, side: %s, model: %s" % (subject, side, model_id))
        if subject == 'rj090242':
            continue

        # step 1. set all the paths
        fmri_dir = os.sep.join((paths[subject]['base'], 
                                paths[subject]['acquisition']))
        write_dir = os.path.join(fmri_dir, result_dir)
        epi_mask = os.path.join(write_dir, 'mask.nii')
        if not os.path.exists(write_dir):
            os.mkdir(write_dir)

        # work on moco data
        sessions = paths[subject]['sessions'] # moco_sessions
        wild_card = paths[subject]['wild_card']
        if side == 'left':
            wild_card = 'r*lh_.gii'
        elif side == 'right':
            wild_card = 'r*rh_.gii'

        # get the images
        fmri_series = [glob.glob(
                os.path.join(fmri_dir, session, wild_card))[0]
                       for session in sessions]
            
        # compute the mask
        if side == False:
            mask_array = compute_mask_files(fmri_series[0], epi_mask, True, 
                                            inf_threshold, sup_threshold)[0]
        # get the contrasts
        n_reg = make_dmtx(frametimes, add_regs=reg_matrix,
                          drift_model=drift_model, hfcut=hfcut).matrix.shape[1]

        contrasts, all_reg = make_contrasts(sessions, n_reg)
        contrast_obj = {}

        for (sess, (session, fmri_data)) in enumerate(zip(
                sessions, fmri_series)):
            # create design matrices
            reg = all_reg[2 * sess: 2 * sess + 2] # fixme
            design_matrix = make_dmtx(
                frametimes, add_regs=reg_matrix, add_reg_names=reg,
                drift_model=drift_model, hfcut=hfcut)

            # plot the design matrix
            ax = design_matrix.show()
            ax.set_position([.05, .25, .9, .65])
            ax.set_title('Design matrix')
            pylab.savefig(os.path.join(write_dir, 'design_matrix_%s.png') %\
                              session)
            # get the data
            if side == False:
                Y, _ = data_scaling(load(fmri_data).get_data()[mask_array].T)
                affine = load(fmri_data).get_affine()
            else:
                Y = np.array([x.data for x in read(fmri_data).darrays])
                Y, _ = data_scaling(Y)
                if sess == 0: # do not redefine the mask later !
                    mask_array = np.var(Y, 0) > 0
                Y = Y[:, mask_array]

            # fit the glm
            print 'Fitting a GLM (this takes time)...'
            result = GeneralLinearModel(design_matrix.matrix)
            result.fit(Y, model='ar1', steps=100)

            for contrast_id, contrast_val in contrasts.iteritems():
                if (contrast_val[session] == 0).all():
                    continue
                contrast_obj[contrast_id] =\
                    result.contrast(contrast_val[session])

        print 'Computing contrasts...'
        # In fact, this only writes the output images/textures
        # define the effects of interest
        F_ = np.zeros(mask_array.sum())
        for index, contrast_id in enumerate(contrasts):
            F_ += (contrast_obj[contrast_id].stat()) ** 2
            print '  Contrast % 2i out of %i: %s' % (
                index + 1, len(contrasts), contrast_id)
            write_array = mask_array.astype(np.float)
            write_array[mask_array] = contrast_obj[contrast_id].effect.ravel()
            if side == False:
                contrast_image = Nifti1Image(write_array, affine)
                save(contrast_image, os.path.join(write_dir, '%s_con.nii' %
                                                  contrast_id))
            else:
                contrast_path = os.path.join(write_dir, '%s_%s_con.gii' %
                                             (side, contrast_id))
                save_texture(contrast_path, write_array)
        # save the effects of interest
        write_array[mask_array] = st.norm.isf(st.f.sf(F_ / 8, 8, 480))
        if side == False:
            contrast_image = Nifti1Image(write_array, affine)
            save(contrast_image, 
                 os.path.join(write_dir, 'effects_of_interest_z_map.nii'))
        else:
            contrast_path = os.path.join(
                write_dir, '%s_effects_of_interest_z_map.gii' % side)
            save_texture(contrast_path, write_array)
        

#--------------------------------------------------------------------
# Retinotopy specific analysis: phase maps
#--------------------------------------------------------------------
    
for subject in subjects:
    print ('Computing phase maps in subject %s' % subject)
    if side != False and subject == 'rj090242':
        continue # missing data for this subject
    for side in sides:
        paths[subject]["contrasts"] = os.path.join(
            paths[subject]['base'], paths[subject]['acquisition'], result_dir)
    
        # offset_wedge and offset_ring are related to the encoding 
        # of the stimulus
        offset_wedge = 0
        # means that the wedge starts at the right horizontal position
        offset_ring = np.pi
        # means that the ring starts at the middle position
        threshold = 3.
        # threshold of the main effects map, in z-value
        angular_analysis.angular_maps(
            side, paths[subject], all_reg, threshold=threshold,
            size_threshold=50, offset_wedge=offset_wedge,
            offset_ring=offset_ring, smooth=0., do_phase_unwrapping=True)
        
