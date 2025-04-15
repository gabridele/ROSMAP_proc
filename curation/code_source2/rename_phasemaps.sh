#!/bin/bash
cd sub-28680644
 
datalad unlock "120411_03_28680644/Obl_2-echo_GRE_phase_map_e1.json"
mv "120411_03_28680644/Obl_2-echo_GRE_phase_map_e1.json" "sub-28680644_ses-#_acq-BNK20090211_Obl_2-echo_GRE_phase_map_e1.json"

datalad unlock "120411_03_28680644/Obl_2-echo_GRE_phase_map_e1.nii.gz"
mv "120411_03_28680644/Obl_2-echo_GRE_phase_map_e1.nii.gz" "sub-28680644_ses-#_acq-BNK20090211_Obl_2-echo_GRE_phase_map_e1.nii.gz"

datalad save -m "Renamed phase maps according to BIDS format"
cd ..
###########################################################################################################