"""
Script calls Freesurfer in order to obtain functional textures:
This projects the fMRI images to the meshes obtained from freesurfer
then smoothes the textures.

Todo
----
- Integrate it to the `preprocessing.py` script
- Do it with nipy
- checkreg

Note
----
First run: export SUBJECTS_DIR=''
"""
import os
import glob
import commands
#import numpy as np
from nibabel import load, save, Nifti1Image
#from nibabel.freesurfer import mghformat
from retino.angular_analysis import load_texture
import config_retino_3T

FWHM = 5.

# ----------------------------------------------------------------
# data paths
# ----------------------------------------------------------------

_, main_dir, subject_info = config_retino_3T.init_config()

# possibly loop on the subjects
for subject in subject_info.keys():

    # the following ones are necessary in order to get the fmri data paths
    fmri_dir = os.path.join(main_dir, subject, 'fmri')
    fs_dir = os.path.join(main_dir, subject, 't1', subject)
    
    # take the mean image that has been coregistered to the anat
    mean_image = glob.glob(os.path.join(fmri_dir, 'mean*.nii'))[0]
    affine = load(mean_image).get_affine()
    
    # take the fMRI series
    sessions = subject_info[subject]['session_keys']
    fmri_paths = [os.path.join(fmri_dir, '%s_series_%s.nii' %
                               (subject, session)) for session in sessions]

    for fmri_session in fmri_paths:
        # create a temporary image with the right affine transform
        fmri_vol = '/tmp/tmp.nii'
        save(Nifti1Image(load(fmri_session).get_data(), affine), fmri_vol)

        # output names
        # the .gii files will be put in the same directory as the input fMRI
        left_fmri_tex = fmri_session[:-4] + '_lh.gii' 
        right_fmri_tex = fmri_session[:-4] + '_rh.gii'
            
        # run freesrufer command

        commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi lh --projfrac 0.5'
            % (fmri_vol, left_fmri_tex, fs_dir))
        plop = commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi rh --projfrac 0.5'
            % (fmri_vol, right_fmri_tex, fs_dir))

# smooth the data
for subject in subject_info.keys(): 
    fmri_dir = os.path.join(main_dir, subject, 'fmri')
    fs_dir = os.path.join(main_dir, subject, 't1', subject)
    for session in subject_info[subject]['session_keys']:
        for hemi in ['left', 'right']:
            wild_card = '%sh.gii' % hemi[0]
            fmri_series = os.path.join(fmri_dir, '%s_series_%s_%s' %
                                       (subject, session, wild_card))
            tex = load_texture(fmri_series)
            smooth_tex = []
            output = fmri_series[:-4] + '_smooth5.gii'
            plop = commands.getoutput(
                'mri_surf2surf --srcsubject %s --cortex --fwhm %s '\
                    '--hemi %sh --sval %s --tval %s --s %s' % (
                        fs_dir, FWHM, hemi[0], fmri_series, output, fs_dir))
  
