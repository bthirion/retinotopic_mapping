"""
This function gets the data on the Neurospin disks
Should be included in the pre-rpocessing script.
"""
import config_retino_7T
#import resample_anat as rsa

import os
import glob
import commands
import shutil
import numpy as np

from nipype.caching import Memory
from nipype.interfaces import spm
import nipype.interfaces.matlab as matlab
from nibabel import load, save, Nifti1Image
from nipy.labs.datasets import VolumeImg
from nipy.labs.datasets.transforms.affine_utils import get_bounds

###############################################################################
# Set the way matlab should be called at Neurospin:
# we have to hit directly the original mathworks shell script
matlab.MatlabCommand.set_default_matlab_cmd(
    "/neurospin/local/matlab/bin/matlab")
matlab.MatlabCommand.set_default_paths('/i2bm/local/spm8/')
T1_TEMPLATE = '/i2bm/local/spm8/templates/T1.nii'

###############################################################################


##############################################################
# Set the paths, dataset-specific stuff
##############################################################

# database path
data_path, subject_info = config_retino_7T.init_config()

# directory where the analysis will take place
main_dir = '/neurospin/tmp/retino/7T/'
if not os.path.exists(main_dir):
    os.mkdir(main_dir)

# parameters
TR = 2.4

for subject in ['eb120536']:#subject_info.keys():
    subject_dict = subject_info[subject]

    # the next part is automatic, possibly neurospin-specific
    subject_dir = os.path.join(main_dir, subject)
    fmri_dir = os.path.join(subject_dir, 'fmri')
    t1_dir = os.path.join(subject_dir, 't1')
    anat_image = os.path.join(t1_dir, '%s_t1.nii' % subject)

    # create directories where to write the image
    for directory in [subject_dir, fmri_dir, t1_dir]:
        if not os.path.exists(directory):
            os.mkdir(directory)

    # now get the ids of the fMRI scans
    fmri_sessions = {}
    for (session_id, session_index) in subject_dict['session_ids'].items():
        if session_id != 't1':
            fmri_sessions[session_id] = os.path.join(
                fmri_dir, '%s_series_%s.nii' % (subject, session_id))
    fmri_series = fmri_sessions.values()
    fmri_series.sort()

    ##############################################################
    # Get the archive and convert the data using MRIcron
    ##############################################################

    for (session_id, session_index) in subject_dict['session_ids'].items():
        # fetch the data
        archive = glob.glob(os.path.join(
            data_path,
            subject_dict['date'], '%s*' % subject_dict['subject_id'],
            '*%s*/' % session_index))[0]

        # copy dicom folder for this session in subject_dir
        dicom_dir = os.path.join(subject_dir, '%s_dicom/' % session_index)
        if os.path.exists(dicom_dir):
            for x in glob.glob(os.path.join(dicom_dir, '*')):
                os.remove(x)
            os.removedirs(dicom_dir)

        shutil.copytree(archive, dicom_dir)
        dcm_files = glob.glob(dicom_dir + '*.dcm')[0]

        output = commands.getoutput(
          'dcm2nii -g N -r N -c N -f N -i Y -d N -p N -e N -o %s %s' % \
            (dicom_dir, dcm_files))

        # move the nifti file to the expected place
        output_idx, = np.where([subject in x for x  in output.split('\n')])
        output_idx = output_idx[-1]
        nifti_file = output.split('\n')[output_idx].split(' ')[-1].split(
                            '>')[-1]
        nifti_file = os.path.join(dicom_dir, nifti_file.split('/')[-1])
        if session_id == 't1':
            print nifti_file, anat_image
            shutil.move(nifti_file, anat_image)
        else:
            print nifti_file, fmri_sessions[session_id]
            shutil.move(nifti_file, fmri_sessions[session_id])
        
        # remove the dicom dirs
        for x in glob.glob(os.path.join(dicom_dir, '*')):
            os.remove(x)
        os.removedirs(dicom_dir)
    
    ##############################################################
    # Preprocessing
    ##############################################################

    ##############################################################
    # Anatomical segmentation (White/Grey matter)
    mem = Memory(base_dir=subject_dir)
    seg = mem.cache(spm.Segment)
    out_seg = seg(data=anat_image,
                  gm_output_type=[True, True, True],
                  wm_output_type=[True, True, True],
                  csf_output_type=[True, True, True])
    sn_file = out_seg.outputs.transformation_mat
    inv_sn_file = out_seg.outputs.inverse_transformation_mat
    gm_image = out_seg.outputs.normalized_gm_image
    native_gm_image = out_seg.outputs.native_gm_image

    shutil.copyfile(native_gm_image, os.path.join(t1_dir,
        '%s_gm_image.nii' % subject))

    ##############################################################
    #  Slice timing correction
    # this can be skipped 
    n_slices = load(fmri_series[0]).get_shape()[2]
    slice_order = range(1, 1 + np.uint8(n_slices))
    mem = Memory(base_dir=main_dir)
    st = mem.cache(spm.SliceTiming)
    st_result = st(in_files=fmri_series,
                   num_slices=n_slices,
                   time_repetition=TR,
                   time_acquisition=TR - TR / n_slices,
                   slice_order=slice_order,
                   ref_slice=n_slices / 2,
                   out_prefix='t')

    time_corrected_files = []
    for stf in st_result.outputs.timecorrected_files:
        time_corrected_files.append(os.path.join(fmri_dir,
                            os.path.basename(stf)))
    
    ##############################################################
    # realign the data just to get the mean image
    realign = mem.cache(spm.Realign)
    realign_result = realign(
        in_files=st_result.outputs.timecorrected_files,
        register_to_mean=True, jobtype='estwrite', out_prefix='r')
    mean_img = realign_result.outputs.mean_image
    shutil.copyfile(mean_img,
                os.path.join(fmri_dir, os.path.basename(mean_img)))

    realigned_files = []
    for rf in realign_result.outputs.realigned_files:
        realigned_files.append(
            os.path.join(fmri_dir, os.path.basename(rf)))
        shutil.copyfile(rf, os.path.join(fmri_dir, os.path.basename(rf)))
    
    for rpf in realign_result.outputs.realignment_parameters:
        shutil.copyfile(rf, os.path.join(fmri_dir, os.path.basename(rpf)))
    
    ##############################################################
    # Coregister the anatomy to the mean functional image 
    #    (only rewrite the header)

    coreg = mem.cache(spm.Coregister)
   
    # get the mean functionnal image to coregister the data on the anat
    mean_image = glob.glob(os.path.join(fmri_dir, 'mean*.nii'))[0]
    
    # does the coregistration of the time corrected+realigned fmri series
    coreg_result = coreg(target=mean_image, source=anat_image,
                         jobtype='estimate')

    """
    ##############################################################
    # Run freesurfer segmentation
    from nipype.interfaces.freesurfer import ReconAll
    reconall =  mem.cache(ReconAll)
    recon_result = reconall(subject_id = subject, 
                            directive='all', 
                            subjects_dir = t1_dir,
                            T1_files = anat_image)
    """
