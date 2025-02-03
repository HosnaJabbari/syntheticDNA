import requests
import os
from concurrent.futures import ThreadPoolExecutor #to dl simultaneously

def download_pdb_or_cif(pdb_id, output_dir, downloaded_log):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb" #URL
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    cif_path = os.path.join(output_dir, f"{pdb_id}.cif")
    try:
        response = requests.get(pdb_url, timeout=10)
        if response.status_code == 200: # if success
            with open(pdb_path, 'wb') as f: #binary cuz dled over http
                f.write(response.content)
            print(f"{pdb_id} PDB dled")
            log_downloaded_file(pdb_id, downloaded_log)  # Log success download
            return pdb_path
        else:
            print(f"Failed {pdb_id} PDB. now CIF...")
        # Try CIF
        response = requests.get(cif_url, timeout=10)
        if response.status_code == 200:
            with open(cif_path, 'wb') as f:
                f.write(response.content)
            print(f"CIF for {pdb_id} DLed")
            log_downloaded_file(pdb_id, downloaded_log)
            with open("cif_downloads.txt", 'a') as log_file: # to log cifs in case necessary in future
                log_file.write(f"{pdb_id}.cif\n")
            return cif_path
        else:
            print(f"Failed {pdb_id} CIF")
            with open("failed_downloads.txt", 'a') as fail_log:
                fail_log.write(f"{pdb_id}\n")
            return None
    except Exception as e:
        print(f"Error {pdb_id}: {e}")
        with open("failed_downloads.txt", 'a') as fail_log:
            fail_log.write(f"{pdb_id}: {e}\n")
        return None

def log_downloaded_file(pdb_id, downloaded_log):#Log success dled PDB/CIF ID to log file
    with open(downloaded_log, 'a') as log_file:
        log_file.write(f"{pdb_id}\n")

def read_pdb_ids(file_path): #read ids from a single line comma-sep format (copied from rcsb output) and return list
    with open(file_path, 'r') as file:
        line = file.readline()
        pdb_ids = line.strip().split(',')  # Split commas & remove whitespace
    return pdb_ids

def read_logged_downloads(downloaded_log):#to read the dl log to avoid already dled files in case the list grows from different databases and to consider concurrent dl process
    if os.path.exists(downloaded_log):
        with open(downloaded_log, 'r') as log_file:
            return set(line.strip() for line in log_file.readlines())
    return set()

def download_all_pdbs(pdb_ids, output_directory, downloaded_log): #multithreading (concurrent)
    os.makedirs(output_directory, exist_ok=True) #ok if its there already:D
    downloaded_files = read_logged_downloads(downloaded_log)
    print(f"Skipping {len(downloaded_files)} already dled.")
    # Filter PDB IDs to download
    pdb_ids_to_download = []
    for pdb_id in pdb_ids:
        if pdb_id not in downloaded_files:
            pdb_ids_to_download.append(pdb_id)
    with ThreadPoolExecutor(max_workers=10) as executor:  #use 10 threads
        futures = {executor.submit(download_pdb_or_cif, pdb_id, output_directory, downloaded_log): pdb_id for pdb_id in pdb_ids_to_download}
        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Error {futures[future]} {e}")

if __name__ == "__main__":
    pdb_ids = read_pdb_ids("all_ids_comma_delimited.txt")
    print(f"Loaded {len(pdb_ids)} PDB IDs")
    output_directory = "RCSB_pdb_files"
    downloaded_log = "downloaded_pdbs.log"
    download_all_pdbs(pdb_ids, output_directory, downloaded_log)
    print("complete")