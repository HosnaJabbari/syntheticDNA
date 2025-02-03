import os
import json

def compare_pdb_ids(ecRBPome_ids, rcsb_ids):
#Compare ecRBPome IDs with RCSB IDs to find common and unique IDs.
    common_ids = ecRBPome_ids & rcsb_ids  # intersection
    unique_to_ecRBPome = ecRBPome_ids - rcsb_ids  # ecRBPome unique
    unique_to_RCSB = rcsb_ids - ecRBPome_ids  # RCSB unique
    return common_ids, unique_to_ecRBPome, unique_to_RCSB

def read_ids_in_file(file_path): # read PDB IDs from file into a set
    with open(file_path, 'r') as file:
        unique_ids = set()
        for line in file:
            ids = line.strip().split(',') #split by comma if there is, otherwise just strip the line
            for id in ids:#add each ID to set
                unique_ids.update([id.strip()])
        return unique_ids

if __name__ == "__main__":
    ecRBPome_file = "ecRBPome_PDB_IDs_unique.txt"
    ecRBPome_ids = read_ids_in_file(ecRBPome_file)

    #read RCSB from config.json
    with open("config.json", 'r') as config_file:
        config_data = json.load(config_file)
    combined_rcsb_file = config_data.get("combined_output_file")

    if not combined_rcsb_file:
        raise ValueError("combined RCSB file not found in config.json")

    rcsb_ids = read_ids_in_file(combined_rcsb_file)
    common_ids, unique_to_ecRBPome, unique_to_RCSB = compare_pdb_ids(ecRBPome_ids, rcsb_ids)
    #logs
    with open("comparison_log_common.txt", 'w') as common_file:
        common_file.write(f"{len(common_ids)} common PDB IDs:\n")
        common_file.write("\n".join(common_ids))
    
    with open("comparison_log_unique_to_ecRBPome.txt", 'w') as unique_ec_file:
        unique_ec_file.write(f"{len(unique_to_ecRBPome)} unique to ecRBPome PDBs,\n")
        unique_ec_file.write("\n".join(unique_to_ecRBPome))
    
    with open("comparison_log_unique_to_RCSB.txt", 'w') as unique_rcsb_file:
        unique_rcsb_file.write(f"{len(unique_to_RCSB)} unique to RCSB PDBs,\n")
        unique_rcsb_file.write("\n".join(unique_to_RCSB))
    with open("all_ids_comma_delimited.txt", 'w') as all_ids_file:
        all_ids = sorted(ecRBPome_ids | rcsb_ids)  # union (no dups),sorted
        all_ids_file.write(",".join(all_ids))
    print("comparison done and saved in log files")
