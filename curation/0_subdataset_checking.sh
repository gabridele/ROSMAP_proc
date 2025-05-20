#!/bin/bash

for sub in sub-*; do
    echo "Processing $sub"

    datalad create -d . -f $sub

    datalad save -d $sub -m 'checked in scans in subject subdataset' 

done