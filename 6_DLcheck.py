import os

def comp_pdb_with_ids(files, pdb_ids):  # compare file names with IDs
    files_in_directory = set(os.path.splitext(file_name)[0] for file_name in files) #set to comp fast
    missing_ids = pdb_ids - files_in_directory  # IDs in file but not in directory
    extra_files = files_in_directory - pdb_ids  # Files in directory but not in file
    return missing_ids, extra_files

with open("all_ids_comma_delimited.txt", 'r') as file: #extracted ids
    pdb_ids = set(file.readline().strip().split(','))

files = os.listdir("RCSB_pdb_files") #dled
missing_ids, extra_files = comp_pdb_with_ids(files, pdb_ids)

if missing_ids:
    print("Missing in directory")
    print(", ".join(missing_ids))
else:
    print("No missing")

if extra_files:
    print("Extra in directory")
    print(", ".join(extra_files))
else:
    print("No extra")