"""
This file prepares some paths for the analysis of the lookloc_data.
It is simply a representation of where the data is on the file system
"""

import os.path as op
import glob

# ----------------------------------------------------------------
# data paths
# ----------------------------------------------------------------



def make_paths():

    db_path = "/volatile/thirion/data/retino_7T" 
    fs_db = None 
    subjects = [ 'gm110134', 'jh100405', 'kr080082', 'vr100551', 'rj090242', 
                 'td110140'] 
    acquisition = ''
    sessions = ['session_1', 'session_2', 'session_3', 'session_4']
    wild_card = "rf*.nii"

    paths = {}

    for subject in subjects:
        paths[subject] = {
            'base': op.join(db_path, subject),
            'acquisition': acquisition,
            'sessions': sessions,
            'mean_image': \
                glob.glob(op.sep.join((
                            db_path, subject, acquisition, sessions[0],
                            'meanf*.img')))[0],
            'epi_mask':\
                glob.glob(op.sep.join((
                            db_path, subject, acquisition, sessions[0],
                            'mask.nii')))[0],
            'motion':\
                glob.glob(op.sep.join((
                            db_path, subject, acquisition, sessions[0],
                            'rp_f*.txt')))[0],
            'wild_card' : wild_card
            }
        if fs_db is None:
            fs_db = db_path
    
        paths[subject]['t1'] = glob.glob('%s/%s/t1/*.img' % 
                                                (fs_db, subject))[0]
        # freesurfer stuff
        paths[subject]['left_mesh'] = '%s/%s/t1/%s/surf/lh.white.gii' % (
            fs_db, subject, subject)
        paths[subject]['right_mesh'] = '%s/%s/t1/%s/surf/rh.white.gii' % (
            fs_db, subject, subject)
        paths[subject]['left_inflated'] = \
            '%s/%s/t1/%s/surf/lh.inflated.gii' % (fs_db, subject, subject)
        paths[subject]['right_inflated'] = \
            '%s/%s/t1/%s/surf/rh.inflated.gii' % (fs_db, subject, subject)
        paths[subject]['left_spherical'] = \
            '%s/%s/t1/%s/surf/lh.sphere.reg.gii' % (fs_db, subject, subject)
        paths[subject]['right_spherical'] =\
            '%s/%s/t1/%s/surf/rh.sphere.reg.gii' % (fs_db, subject, subject)
        paths[subject]['reg_file'] = '%s/%s/t1/%s/mri/orig.mgz' % (
            fs_db, subject, subject)
        paths[subject]['fs_subj'] = subject
        

    return paths


def init_config():
    """Initialize texture decoding experiement specific variables

    Returns
    =======
    main_dir: string, directory containing texture decoding experiment mri data
    subjects_information: 
        dictionary, subject informations : subject id, session ids,
        session date, stimuli set (1 or 2)
    """
    # data repository
    data_path = '/neurospin/acquisition/database/Investigational_Device_7T/'

    subject_info = {
        'eb120536': {
            'folder': 'eb120536-934_001',
            'subject_id': 'eb120536',
            'session_ids': {
                't1': '000017_t1-mpr-tra-iso1.0mm',
                'ring_neg': '000014_cmrr-mbep2d-iso1mm-p2-MB2-80sl-132rep',
                'ring_pos': '000020_cmrr-mbep2d-iso1mm-p2-MB2-80sl-132rep',
                'wedge_pos': '000015_cmrr-mbep2d-iso1mm-p2-MB2-80sl-132rep',
                'wedge_neg': '000016_cmrr-mbep2d-iso1mm-p2-MB2-80sl-132rep',
                'wedge_pos2': '000021_cmrr-mbep2d-iso1mm-p2-MB2-80sl-132rep'},
            'date': '20121128',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'gm110134':{
            'folder': 'gm110134-778_001',
            'subject_id': 'gm110134',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'ring_pos': '000006_MoCoSeries',
                'ring_neg': '000008_MoCoSeries',
                'wedge_pos': '000010_MoCoSeries',
                'wedge_neg': '000012_MoCoSeries'},
            'date': '20110914',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'jh100405':{
            'folder': 'jh100405-779_001',
            'subject_id': 'jh100405',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'ring_neg': '000012_MoCoSeries',
                'ring_pos': '000014_MoCoSeries',
                'wedge_pos': '000017_MoCoSeries',
                'wedge_neg': '000019_MoCoSeries'},
            'date': '20110914',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'vr100551':{
            'folder': 'vr100551',
            'subject_id': 'vr100551-786_001',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'ring_neg': '000009_MoCoSeries',
                'ring_pos': '000011_MoCoSeries',
                'wedge_pos': '000013_MoCoSeries',
                'wedge_neg': '000015_MoCoSeries'},
            'date': '20110928',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'kr080082':{
            'folder': 'kr080082-787_001',
            'subject_id': 'kr080082',
            'session_ids': {
                't1': '000003_mprage-sag-T1-160sl',
                'ring_neg': '000007_MoCoSeries',
                'ring_pos': '000009_MoCoSeries',
                'wedge_pos': '000015_MoCoSeries',
                'wedge_neg': '000017_MoCoSeries'},
            'date': '20110928',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'rj090242':{
            'folder': 'rj090242-790_001',
            'subject_id': 'rj090242',
            'session_ids': {
                't1': '000003_mprage-sag-T1-160sl',
                'ring_neg': '000007_MoCoSeries',
                'ring_pos': '000009_MoCoSeries',
                'wedge_pos': '000011_MoCoSeries',
                'wedge_neg': '000013_MoCoSeries'},
            'date': '20111012',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        'td110140':{
            'folder': 'td110140-791_001',
            'subject_id': 'td110140',
            'session_ids': {
                't1': '000003_mprage-sag-T1-160sl',
                'ring_neg': '000011_MoCoSeries',
                'ring_pos': '000013_MoCoSeries',
                'wedge_pos': '000015_MoCoSeries',
                'wedge_neg': '000017_MoCoSeries'},
            # CAVEAT: there is also the series 007, not sure of what happened
            'date': '20111012',
            'protocol': 'ring + wedge',
            'scanner': '7T'},
        }    
    return data_path, subject_info

if __name__ == 'main':
    make_paths()
