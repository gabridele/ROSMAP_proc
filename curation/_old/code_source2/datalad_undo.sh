#!/bin/bash

function datalad_undo {
    dir="$1"
    echo "Processing directory: $dir"
    cd ./$dir

    # Get the latest commit message
    latest_commit_message=$(git log -1 --pretty=%B)

    # Check if the message contains '[...]'
    if [[ "$latest_commit_message" == *"3D_MPRAGE"* ]]; then

        search_term="populated subject"
        sha=$(git log --pretty=oneline | grep "$search_term" | awk '{print $1}')
        position=$(git rev-list --count $sha..HEAD)
        echo "Resetting to HEAD~$position"
        datalad unlock ./*
        git reset --hard HEAD~$position

    else
        echo "Commit message in $dir doesnt contain $latest_commit_message."
        
    fi
    cd ..

}


N=1
(
while IFS= read -r dir; do
    ((i=i%N)); ((i++==0)) && wait
    datalad_undo "${dir%%,*}" &
done < "code/subs_batch2.txt"
wait
)
