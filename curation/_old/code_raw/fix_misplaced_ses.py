import os
import datalad.api as dl

subjects = [
    "sub-15938020",
    "sub-18455382",
    "sub-24644776",
    "sub-27083391",
    "sub-44842532",
    "sub-48480640",
    "sub-49094578",
    "sub-50404037",
    "sub-52052115",
    "sub-65778949",
    "sub-65925304",
    "sub-65952361",
    "sub-67946467",
    "sub-73177635",
    "sub-78010047",
    "sub-79692387",
    "sub-81185560",
    "sub-81810992",
    "sub-83216408",
    "sub-84896566",
    "sub-87585034",
    "sub-90805280",
]

for subject in subjects:
    dl.unlock(subject)

def rename_dir(source, target):
    """
    Rename a ses- directory and the ses- string in filename accordingly.
    Args:
        source (str): Source directory to rename.
        target (str): Target directory name.
    """
    try:
        os.rename(source, target)
        print(f"Renamed {source} to {target}")
        
        # Update file names inside the directory if the source string matches
        for root, dirs, files in os.walk(target):
            for dir_name in dirs:
                dir_path = os.path.join(root, dir_name)
                for sub_root, sub_dirs, sub_files in os.walk(dir_path):
                    for file in sub_files:
                        if source in file:
                            old_file_path = os.path.join(sub_root, file)
                            new_file_name = file.replace(source, target)
                            new_file_path = os.path.join(sub_root, new_file_name)
                            os.rename(old_file_path, new_file_path)
                            print(f"Renamed {old_file_path} to {new_file_path}")
    except OSError as e:
        print(f"Error renaming {source} to {target}: {e}")

def rename_ses(dir, match, target):
    """
    Rename ses- sring in filename files, according to specified pattern.
    Args:
        dir (str): Directory containing the files.
        match (str): String to match in the file names.
        target (str): String to replace the matched string.
    """
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)  # Full path to the file
        if match in file:
            new_file = file.replace(match, target)
            new_file_path = os.path.join(dir, new_file)  # Full path for the new file
            os.rename(file_path, new_file_path)  # Use full paths
            print(f"Renamed {file_path} to {new_file_path}")
        else:
            print(f"No match found for {file_path}")

def rm_duplicates(dir, match):
    """
    Remove duplicate files that contain the specified match in their names.
    
    Args:
        dir (str): Directory containing the files.
        match (str): String to match in the file names.
    """
    for sub_dir in os.listdir(dir):
        sub_dir_path = os.path.join(dir, sub_dir)
        if os.path.isdir(sub_dir_path):
            for file in os.listdir(sub_dir_path):
                file_path = os.path.join(sub_dir_path, file)
                print(f"Checking file: {file_path} against match: {match}")  # Debugging output
                if match in file:
                    os.remove(file_path)
                    print(f"Removed duplicate file: {file_path}")
                else:
                    print(f"No match found for {file_path}")

def move(source, target, match):
    """
    Move all files matching the pattern from source to target.

    Args:
        source (str): Source directory containing the files.
        target (str): Target directory to move the files to.
        match (str): String to match in the file names.
    """
    if not os.path.exists(target):
        os.makedirs(target)
        print(f"Created directory: {target}")
    for file in os.listdir(source):
        if match in file:
            src_path = os.path.join(source, file)
            dest_path = os.path.join(target, file)
            os.rename(src_path, dest_path)
            print(f"Moved {file} from {source} to {target}")
        else:
            print(f"No match found for {file}")

def add_BNK(dir):
    """
    Add 'BNK' site to the acq- string in file name.

    Args:
        dir (str): Directory containing the files.
    """
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)  # Full path to the file
        if "acq-2009" in file:
            new_file = file.replace("acq-2009", "acq-BNK2009")
            new_file_path = os.path.join(dir, new_file)  # Full path for the new file
            os.rename(file_path, new_file_path)  # Use full paths
            print(f"Renamed {file_path} to {new_file_path}")
        else:
            print(f"No match found for {file_path}")

def rename_acq(dir, match, target):
    """
    Rename acq string in filenames according to the specified pattern.

    Args:
        dir (str): Directory containing the files.
        match (str): String to match in the file names.
        target (str): String to replace the matched string.
    """
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)  # Full path to the file
        if match in file:
            new_file = file.replace(match, target)
            new_file_path = os.path.join(dir, new_file)  # Full path for the new file
            os.rename(file_path, new_file_path)  # Use full paths
            print(f"Renamed {file_path} to {new_file_path}")
        else:
            print(f"No match found for {file_path}")

def main():
    """
    Main function to handle renaming and moving files.
    """

    ## current directory is 'raw', which contains all the sub-*
    rename_dir("sub-48480640/ses-2", "sub-48480640/ses-3")

    rename_dir("sub-65952361/ses-2", "sub-65952361/ses-3")
    rename_dir("sub-65952361/ses-1", "sub-65952361/ses-2")
    
    rename_dir("sub-67946467/ses-2", "sub-67946467/ses-3")

    rename_dir("sub-78010047/ses-3", "sub-78010047/ses-4")
    rename_dir("sub-78010047/ses-2", "sub-78010047/ses-3")
    rename_dir("sub-78010047/ses-1", "sub-78010047/ses-2")
    rename_dir("sub-78010047/ses-0", "sub-78010047/ses-1")
    
    rename_dir("sub-79692387/ses-2", "sub-79692387/ses-3")

    rename_dir("sub-81185560/ses-4", "sub-81185560/ses-5")
    rename_dir("sub-81185560/ses-3", "sub-81185560/ses-4")
    rename_dir("sub-81185560/ses-2", "sub-81185560/ses-3")

    rename_dir("sub-81810992/ses-3", "sub-81810992/ses-4")
    rename_dir("sub-81810992/ses-2", "sub-81810992/ses-3")
    rename_dir("sub-81810992/ses-1", "sub-81810992/ses-2")
    rename_dir("sub-81810992/ses-0", "sub-81810992/ses-1")
    
    rename_dir("sub-83216408/ses-3", "sub-83216408/ses-4")
    rename_dir("sub-83216408/ses-2", "sub-83216408/ses-3")
    rename_dir("sub-83216408/ses-1", "sub-83216408/ses-2")
    rename_dir("sub-83216408/ses-0", "sub-83216408/ses-1")

    rename_dir("sub-87585034/ses-2", "sub-87585034/ses-3")
## rename content inside too


    modalities = ["anat", "fmap", "func"]
    for modality in modalities:
        
        move(f"sub-15938020/ses-0/{modality}", f"sub-15938020/ses-1/{modality}", "MG")
        rename_ses(f"sub-15938020/ses-1/{modality}", "ses-0", "ses-1")
        
        
        add_BNK(f"sub-18455382/ses-0/{modality}")
        
        move(f"sub-18455382/ses-1/{modality}", f"sub-18455382/ses-2/{modality}", "MG")
        rename_ses(f"sub-18455382/ses-2/{modality}", "ses-1", "ses-2")


        add_BNK(f"sub-24644776/ses-0/{modality}")

        move(f"sub-24644776/ses-1/{modality}", f"sub-24644776/ses-2/{modality}", "MG")
        rename_ses(f"sub-24644776/ses-2/{modality}", "ses-1", "ses-2")

        rm_duplicates(f"sub-27083391/ses-1/{modality}", "acq-2009")
        ## remove current ses-0
        ## move bnks in ses-1 to ses-0 and rename the files accordingly

        move(f"sub-44842532/ses-1/{modality}", f"sub-44842532/ses-2/{modality}", "MG")
        rename_ses(f"sub-44842532/ses-2/{modality}", "ses-1", "ses-2")


        rename_acq(f"sub-48480640/ses-3/{modality}", "acq-20120501", "acq-MG20150715")

        rename_acq(f"sub-48480640/ses-3/{modality}", "acq-MG20120501", "acq-MG20150715")

        move(f"sub-48480640/ses-1/{modality}", f"sub-48480640/ses-2/{modality}", "MG")
        rename_ses(f"sub-48480640/ses-2/{modality}", "ses-1", "ses-2")

        move(f"sub-49094578/ses-1/{modality}", f"sub-49094578/ses-2/{modality}", "MG")
        rename_ses(f"sub-49094578/ses-2/{modality}", "ses-1", "ses-2")

        rm_duplicates(f"sub-50404037/ses-0/{modality}", "acq-2009")

        move(f"sub-50404037/ses-0/{modality}", f"sub-50404037/ses-1/{modality}", "MG")
        rename_ses(f"sub-50404037/ses-1/{modality}", "ses-0", "ses-1")


        move(f"sub-52052115/ses-0/{modality}", f"sub-52052115/ses-1/{modality}", "MG")
        rename_ses(f"sub-52052115/ses-1/{modality}", "ses-0", "ses-1")


        rm_duplicates(f"sub-65778949/ses-1/{modality}", "acq-2009")

        move(f"sub-65778949/ses-1/{modality}", f"sub-65778949/ses-2/{modality}", "MG")
        rename_ses(f"sub-65778949/ses-2/{modality}", "ses-1", "ses-2")

        add_BNK(f"sub-65925304/ses-0/{modality}")
        rm_duplicates(f"sub-65925304/ses-1/{modality}", "acq-2009")
        
        rename_acq(f"sub-65925304/ses-2/{modality}", "ses-2_task-rest_acq-MG20120501", "ses-2_task-rest_acq-acq-MG20160627")
        move(f"sub-65925304/ses-2/{modality}", f"sub-65925304/ses-3/{modality}", "MG")
        rename_ses(f"sub-65925304/ses-3/{modality}", "ses-2", "ses-3")
        
        move(f"sub-65925304/ses-1/{modality}", f"sub-65925304/ses-2/{modality}", "MG")
        rename_ses(f"sub-65925304/ses-2/{modality}", "ses-1", "ses-2")

        add_BNK(f"sub-65952361/ses-0/{modality}")

        move(f"sub-65952361/ses-0/{modality}", f"sub-65952361/ses-1/{modality}", "MG")
        rename_acq(f"sub-65952361/ses-3/{modality}", "acq-MG20120501", "acq-MG20160627")
        rename_acq(f"sub-65952361/ses-3/{modality}", "acq-20120501", "acq-MG20160627")
        
        
        add_BNK(f"sub-67946467/ses-0/{modality}")
        rm_duplicates(f"sub-67946467/ses-1/{modality}", "acq-2009")

        move(f"sub-67946467/ses-1/{modality}", f"sub-67946467/ses-2/{modality}", "MG")
        rename_ses(f"sub-67946467/ses-2/{modality}", "ses-1", "ses-2")


        move(f"sub-73177635/ses-0/{modality}", f"sub-73177635/ses-2/{modality}", "MG")
        rename_ses(f"sub-73177635/ses-2/{modality}", "ses-0", "ses-2")


        move(f"sub-78010047/ses-1/{modality}", f"sub-78010047/ses-0/{modality}", "BNK")
        rename_acq(f"sub-78010047/ses-2/{modality}", "acq-20120501", "acq-MG20120501")


        rm_duplicates(f"sub-79692387/ses-1/{modality}", "acq-2009")
        rename_acq(f"sub-79692387/ses-3/{modality}", "acq-20120501", "acq-MG20160627")
        rename_acq(f"sub-79692387/ses-3/{modality}", "acq-MG20120501", "acq-MG20160627")

        move(f"sub-79692387/ses-1/{modality}", f"sub-79692387/ses-2/{modality}", "MG")


        rm_duplicates(f"sub-81185560/ses-1/{modality}", "acq-2009")

        move(f"sub-81185560/ses-1/{modality}", f"sub-81185560/ses-2/{modality}", "MG")
        rename_ses(f"sub-81185560/ses-2/{modality}", "ses-1", "ses-2")


        rm_duplicates(f"sub-81810992/ses-1/{modality}", "acq-2009")

        rename_acq(f"sub-81810992/ses-3/{modality}", "acq-20120501", "acq-MG20160627")
        rename_acq(f"sub-81810992/ses-3/{modality}", "acq-MG20120501", "acq-MG20160627")

        move(f"sub-81810992/ses-1/{modality}", f"sub-81810992/ses-0/{modality}", "BNK")
        rename_ses(f"sub-81810992/ses-0/{modality}", "ses-1", "ses-0")

        move(f"sub-83216408/ses-1/{modality}", f"sub-83216408/ses-0/{modality}", "BNK")
        rename_ses(f"sub-83216408/ses-0/{modality}", "ses-1", "ses-0")
        
        rename_acq(f"sub-83216408/ses-2/{modality}", "acq-20120501", "acq-MG20150715")
        rename_acq(f"sub-83216408/ses-2/{modality}", "acq-MG20120501", "acq-MG20150715")

        rename_acq(f"sub-83216408/ses-3/{modality}", "acq-20150715", "acq-MG20160627")
        rename_acq(f"sub-83216408/ses-3/{modality}", "acq-MG20150715", "acq-MG20160627")

        rm_duplicates(f"sub-84896566/ses-0/{modality}", "acq-2009")
        move(f"sub-84896566/ses-0/{modality}", f"sub-84896566/ses-1/{modality}", "MG")
        rename_ses(f"sub-84896566/ses-1/{modality}", "ses-0", "ses-1")


        rm_duplicates(f"sub-87585034/ses-1/{modality}", "acq-2009")

        rename_acq(f"sub-87585034/ses-3/{modality}", "acq-20120501", "acq-MG20160627")
        rename_acq(f"sub-87585034/ses-3/{modality}", "acq-MG20120501", "acq-MG20160627")

        move(f"sub-87585034/ses-1/{modality}", f"sub-87585034/ses-2/{modality}", "MG")


        rm_duplicates(f"sub-90805280/ses-0/{modality}", "acq-2009")
        move(f"sub-90805280/ses-0/{modality}", f"sub-90805280/ses-2/{modality}", "MG")

if __name__ == "__main__":
    main()
    print("All operations completed.")
