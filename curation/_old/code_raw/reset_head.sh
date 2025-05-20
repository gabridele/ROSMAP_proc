#!/bin/bash

CWD=$(pwd)

for sub in "$CWD"/sub-*; do
    cd $sub
    git reset --hard HEAD
    cd ..
done
