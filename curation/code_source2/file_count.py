import os

def count_files_in_subfolders(base_path):
    """Counts the number of files in each subfolder and generates a report.
    Args:
        base_path (str): The base directory containing subfolders.
    """
    report_file = "code/report.txt"
    with open(report_file, "w") as report:
        for folder in os.listdir(base_path):
            folder_path = os.path.join(base_path, folder)
            if os.path.isdir(folder_path):
                for subfolder in os.listdir(folder_path):
                    subfolder_path = os.path.join(folder_path, subfolder)
                    if os.path.isdir(subfolder_path):
                        if subfolder.startswith("."):
                            continue
                        file_count = len([f for f in os.listdir(subfolder_path) if os.path.isfile(os.path.join(subfolder_path, f))])
                        if file_count != 10:
                            report.write(f"{subfolder_path} has {file_count} files (Expected: 10)\n")

if __name__ == "__main__":
    base_directory = "."  # Change this to your actual base directory
    count_files_in_subfolders(base_directory)
    print("Check report.txt for results.")
