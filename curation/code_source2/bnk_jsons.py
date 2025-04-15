import os

import json
import sys
import shutil
import numpy as np
import nibabel as nib
sys.path.insert(0,'./')
import datalad.api as dl
from glob import glob
import pandas
from pathlib import Path

base_json = {'Modality': 'MR',
             'MagneticFieldStrength': 1.5,
             'Manufacturer': 'GE',
             'ManufacturersModelName': 'SIGNA_EXCITE',
             'InstitutionName': 'Bannockburn_Radiology_Center',
             'DeviceSerialNumber': '00000847GURNEEMR',
             'StationName': 'BRCMR',
             'PatientPosition': 'HFS',
             'SoftwareVersions': '11',
             'TaskName':'rest',
             'Instructions':'Eyes closed',
             'EchoTime': 0.0033,
             'RepetitionTime': 2,
             'FlipAngle': 85,
             'PulseSequenceType': '2D Spiral GRE Resting State fMRI',
             'SliceThickness':5,
            'SliceEncodingDirection':'k',
            'SequenceName':'sprlio',
            'SliceTiming':[0.0,0.08,0.16,0.24,0.32,0.4,0.48,0.56,0.64,0.72,0.8,0.88,0.96,1.04,1.12,1.2,1.28,1.36,1.44,1.52,1.6,1.68,1.76,1.84,1.92,2.0]
            }


def write_bnk_json(pth,base_json):
    with open(pth, 'w') as fp:
        json.dump(base_json, fp,sort_keys=True, indent=4)

fmri_files = os.path.join(os.getcwd(), 'sub*', '*', 'sub-*_ses-#_task-rest_acq-BNK20090211_EPI_bold.nii.gz')

single_files = sorted(glob(fmri_files))

for ii, fmri in enumerate(single_files):
    file_path = Path(fmri)
    directory_path = file_path.parent
    sub_id = file_path.parents[1].name

    jpth = fmri.split('.')[0]+'.json'
    write_bnk_json(jpth,base_json)

    dl.save(dataset=sub_id, message=f"Added json sidecar to rs-fmri in {directory_path}")