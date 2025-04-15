import os
import csv

def get_unique_subdirs(directories):
    """
    Get a list of unique subdirectories starting with 'sub-' from the given directories.
    
    Args:
        directories (list): List of directory paths to search for subdirectories.
    
    Returns:
        list: A list of unique subdirectory names starting with 'sub-'.
    """
    unique_subdirs = set()  # Use a set to ensure uniqueness
    
    for directory in directories:
        # Check if the directory exists
        if os.path.exists(directory):
            # Find all subdirectories starting with 'sub-' in the current directory
            subdirs = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d)) and d.startswith("sub-")]
            unique_subdirs.update(subdirs)  # Add the subdirectories to the set
            # Uncomment the line below for debugging purposes
            # print(f"Found {len(subdirs)} subdirectories in {directory}.")
        else:
            # Print a warning if the directory does not exist
            print(f"Warning: Directory {directory} does not exist.")

    return list(unique_subdirs)  # Convert the set back to a list

if __name__ == "__main__":
    # List of directories to search for subdirectories
    directories = ["."]
    
    # Output file to save the list of participants
    subj_list = "participants.tsv" 
    
    # Get the unique subdirectories starting with 'sub-'
    unique_subdirs = get_unique_subdirs(directories)

    # Write the unique subdirectories to a TSV file
    with open(subj_list, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")  # Use tab as the delimiter for TSV
        writer.writerow(["participant_id"])  # Add header row
        for subdir in sorted(unique_subdirs):  # Sort the subdirectories alphabetically
            writer.writerow([subdir])  # Write each subdirectory to the file

    # Print a message indicating the output file
    print(f"Unique subdirectories saved to {subj_list}.")