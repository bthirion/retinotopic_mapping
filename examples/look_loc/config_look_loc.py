"""
This file prepares some paths for the analysis of the lookloc_data.
It is simply a representation of where the data is on the file system
"""

import os.path as op
import glob

# ----------------------------------------------------------------
# data paths
# ----------------------------------------------------------------

db_path = "/volatile/thirion/data/lookloc" # "/volatile/thirion/lookloc" # 
fs_db = None # "/volatile/thirion/fs_db"
subjects = ["sujet04", "sujet10", "sujet11", "sujet12", "sujet13",
            "sujet16", "sujet18", "sujet20"]
wild_card = "rfLOOK*.img"

paths = {}

# freesurfer path
fs_id = {'sujet04': 'sujet04',
         'sujet10': 'sujet10',
         'sujet11': 'sujet11',
         'sujet12': 'sujet12',
         'sujet13': 'sujet13',
         'sujet16': 'sujet16',
         'sujet18': 'sujet18',
         'sujet20': 'sujet20'} 

def make_paths(freesurfer=True, fs_db=fs_db):
    for subject in subjects:
        if subject == 'sujet04':
            acquisition = "S04_2010_04_17_mri14"
            sessions = ["S06_EP_iso2_ring", "S07_EP_iso2_wedge"]
        elif subject == 'sujet10':
            acquisition = "S10_2010_04_11_mri07"
            sessions = ["S03_EP_iso2_ring", "S04_EP_iso2_wedge"]
        elif subject == 'sujet11':
            acquisition = 'S11_2010_04_27_mri15'
            sessions = ['S04_EP_iso2_ring', 'S05_EP_iso2_wedge']
        elif subject == 'sujet12':
            acquisition = 'S12_2010_04_11_mri09'
            sessions = ['S03_EP_iso2_ring', 'S04_EP_iso2_wedge']
        elif subject == 'sujet13':
            acquisition = 'S13_2010_03_31_mri05'
            sessions = ['S03_EP_iso2_ring', 'S04_EP_iso2_wedge']
        elif subject == 'sujet16':
            acquisition = 'S16_2010_06_28_mri18'
            sessions = ['S03_EP_iso2_ring', 'S04_EP_iso2_wedge']
        elif subject == 'sujet18':
            acquisition = 'S18_2010_04_10_mri06'
            sessions = ['S03_EP_iso2_ring', 'S04_EP_iso2_wedge']
        elif subject == 'sujet20':
            acquisition = 'S20_2010_07_04_mri20'
            sessions = ['S03_EP_iso2_ring', 'S04_EP_iso2_wedge']

        paths[subject] = {
        'base': op.join(db_path, subject),
        'acquisition': acquisition,
        'sessions': sessions,
        'mean_image': \
        glob.glob(op.sep.join((db_path, subject, acquisition, sessions[0],
                     'meanf*.img')))[0],
        'wild_card' : wild_card
        }
        if fs_db is None:
            fs_db = db_path
        if freesurfer:
            #freesurfer stuff
            paths[subject]['left_mesh'] = '%s/%s/surf/lh.white.gii' % (
                fs_db, fs_id[subject])
            paths[subject]['right_mesh'] = '%s/%s/surf/rh.white.gii' % (
                fs_db, fs_id[subject])
            paths[subject]['left_inflated'] = '%s/%s/surf/lh.inflated.gii' % (
                fs_db, fs_id[subject])
            paths[subject]['right_inflated'] = '%s/%s/surf/rh.inflated.gii' % (
                fs_db, fs_id[subject])
            paths[subject]['left_spherical'] = '%s/%s/surf/lh.sphere.reg.gii'\
                                               % (fs_db, fs_id[subject])
            paths[subject]['right_spherical'] = '%s/%s/surf/rh.sphere.reg.gii'\
                                                % (fs_db, fs_id[subject])
            paths[subject]['fs_subj'] = fs_id[subject]
        else:
            anat_path = op.sep.join(( db_path, subject,
                        'brainvisa/LOOKLOC/S04/t1mri/default_acquisition/')),
            mesh_path = op.join(anat_path,
                                'default_analysis/segmentation/mesh/'),
            # mesh and kernel paths in brainvisa
            paths[subject]['anat'] = op.join(anat_path, 'S04.nii'),
            paths[subject]['left_mesh'] = op.join(mesh_path, 'S04_Lwhite.gii')
            paths[subject]['right_mesh'] = op.join(mesh_path, 'S04_Rwhite.gii')
            paths[subject]['left_kernel'] = op.join(
                mesh_path, '%s_lh_kernel.ima' % subject)
            paths[subject]['right_kernel'] = op.join(
                mesh_path, '%s_rh_kernel.ima' % subject)

    
    return paths


if __name__ == 'main':
    make_paths()
