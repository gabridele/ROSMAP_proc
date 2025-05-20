#!/bin/bash

# Get the current working directory
CWD=$(pwd)

# Loop through all directories starting with "sub-" in the current working directory
for sub in "$CWD"/sub-*; do
    # Change into the subject directory
    cd $sub

    # Define the search term to look for in the latest commit message
    search_term="Rename acquisition block in BIDS filenames"

    # Get the latest commit message
    latest_commit_message=$(git log -1 --pretty=%B)

    # Check if the latest commit message contains the search term
    if [[ "$latest_commit_message" == *"$search_term"* ]]; then
        echo "Latest commit matches search term. Resetting to HEAD~1"
        # Reset the repository to the previous commit
        git reset --hard HEAD~1
    else
        echo "Latest commit does not match search term. No reset performed."
    fi

    # Return to the parent directory
    cd ..
done
