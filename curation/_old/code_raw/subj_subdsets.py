import os
import datalad.api as dl

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

def create_datalad_datasets(superdataset_path, subdirs):
    """
    Create Datalad datasets for each unique subdirectory.

    Args:
        superdataset_path (str): Path to the superdataset where subdatasets will be created.
        subdirs (list): List of unique subdirectory names to create datasets for.
    """
    # Change to the superdataset directory
    os.chdir(superdataset_path)
    
    for subdir in subdirs:
        subdir_path = os.path.join(superdataset_path, subdir)  # Full path to the subdirectory
        if not os.path.exists(subdir_path):
            # Create a new Datalad dataset for the subdirectory
            dl.create(dataset=superdataset_path, path=subdir)
        else:
            # Skip if the dataset already exists
            print(f"Dataset {subdir} already exists, skipping.")
    
    print("Datalad datasets created successfully.")

if __name__ == "__main__":
    # List of directories to search for subdirectories
    directories = ["../sourcedata/batch_1", "../sourcedata/batch_2", "../sourcedata/batch_3"]
    
    # Path to the superdataset where subdatasets will be created
    superdataset_path = "." 
    
    # Get the unique subdirectories starting with 'sub-'
    unique_subdirs = get_unique_subdirs(directories)
    
    # Create Datalad datasets for the unique subdirectories
    create_datalad_datasets(superdataset_path, unique_subdirs)