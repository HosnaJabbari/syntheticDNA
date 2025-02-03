def remove_duplicates(input_file, output_file):
# Remove duplicates and save to new file.
    unique_ids = set()   # a set to store them
    with open(input_file, 'r') as file:#read the file and store unique lines
        for line in file:
            unique_ids.add(line.strip()) #strip whitespace and add to the set

    with open(output_file, 'w') as file:
        for unique_id in sorted(unique_ids):  # Sort output
            file.write(f"{unique_id}\n")

    print(f"removed duplicates and {len(unique_ids)} IDs saved to {output_file}.")

# input and output
input_file = "ecRBPome_PDB_IDs.txt"
output_file = "ecRBPome_PDB_IDs_unique.txt"
remove_duplicates(input_file, output_file)
