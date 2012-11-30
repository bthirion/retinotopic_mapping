"""
This file prepares some paths for the analysis of the lookloc_data.
It is simply a representation of where the data is on the file system
"""

import os.path as op
import glob

# ----------------------------------------------------------------
# data paths
# ----------------------------------------------------------------


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
    data_path = '/neurospin/acquisition/database/TrioTim'

    # word_decoding experiment:
    # 'ap100009', 'kr080082', 'mr080072', 'vl100318'
    # 'bt080165' did not do it. What about 'ol120056' ?
    # texture_decoding experiment:
    # 
    # old trials:

    subject_info = {
        'kr080082': {
            'folder': 'kr080082-2041_001',
            'subject_id': 'kr080082',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'wedge_pos': '000018_MoCoSeries',
                'wedge_neg': '000020_MoCoSeries'},
            'date': '20100721',
            'protocol': 'wedge',
            'scanner': '3T'},

        'mr080072': {
            'folder': 'mr080072-2037_001',
            'subject_id': 'mr080072',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'wedge_pos': '000020_MoCoSeries',
                'wedge_neg': '000024_MoCoSeries'},
            'date': '20100720',
            'protocol': 'wedge',
            'scanner': '3T'},       

        'vl100318': {
            'folder': 'vl100318-2038_001',
            'subject_id': 'vl100318',
            'session_ids': {
                't1': '000002_mprage-sag-T1-160sl',
                'wedge_pos': '000018_MoCoSeries',
                'wedge_neg': '000020_MoCoSeries'},
            'date': '20100720',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ap100009': {
            'folder': 'ap100009-1789_001',
            'subject_id': 'ap100009',
            'session_ids': {
                't1': '000002_T1-MPRage-Sag',
                'wedge_pos': '000018_MoCoSeries',
                'wedge_neg': '000020_MoCoSeries'},
            'date': '20100505',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ap100009_2': {
            'folder': 'ap100009-3085_001',
            'subject_id': 'ap100009',
            'session_ids': {
                't1': '000020',
                'wedge_pos': '000017',
                'wedge_neg': '000019'},
            'date': '20120509',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ib100049': {
            'folder': '',
            'subject_id': 'ib100049',
            'session_ids': {
                't1': '000003',
                'wedge_pos': '000019',
                'wedge_neg': '000021'},
            'date': '20120515',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ns110383_1': {
            'folder': '',
            'subject_id': 'ns110383',
            'session_ids': {
                't1': '000002',
                'wedge_pos': '000018',
                'wedge_neg': '000020'},
            'date': '20120705',
            'protocol': 'wedge',
            'scanner': '3T'},

        'rm080030': {
            'folder': '',
            'subject_id': 'rm080030',
            'session_ids': {
                't1': '000002',
                'wedge_pos': '000020',
                'wedge_neg': '000022'},
            'date': '20120515',
            'protocol': 'wedge',
            'scanner': '3T'},

        'pf120155_1': {
            'folder': '',
            'subject_id': 'pf120155',
            'session_ids': {
                't1': '000002',
                'wedge_pos': '000018',
                'wedge_neg': '000020'},
            'date': '20120712',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ap100009_3': {
            'folder': '',
            'subject_id': 'ap100009',
            'session_ids': {
                't1': '000002',
                'wedge_pos': '000018',
                'wedge_neg': '000020'},
            'date': '20120822',
            'protocol': 'wedge',
            'scanner': '3T'},

        'ns110383_2': {
            'folder': '',
            'subject_id': 'ns110383',
            'session_ids': {
                't1': '000002',
                'wedge_pos': '000018',
                'wedge_neg': '000020'},
            'date': '20120906',
            'protocol': 'wedge',
            'scanner': '3T'},

        'pf120155_2': {
            'folder': '',
            'subject_id': 'pf120155',
            'session_ids': {
                't1': '000003',
                'wedge_pos': '000019',
                'wedge_neg': '000021'},
            'date': '20120913',
            'protocol': 'wedge',
            'scanner': '3T'},

        }    
    return data_path, subject_info
