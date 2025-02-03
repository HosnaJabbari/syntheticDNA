# xna_project

**1_RCSB_PDBID_Retrieval.py**

  It retrieves PDB IDs for RNA- and DNA-binding proteins from RCSB PDB. It performs a structured query using the RCSB API and filters results to include only structures from E.Coli as source organism.
  
  The script searches for two terms of "RNA Binding Protein" and "DNA Binding Protein" extracts the corresponding PDB IDs and saves them into separate text files (RCSB_RNA_Binding_Protein.txt and RCSB_DNA_Binding_Protein.txt). Moreover, all retrieved PDB IDs are combined and stored in a single file (RCSB_combined_PDB_IDs.txt) to make downstream processing easier. The script also generates a configuration file (config.json) which lists the output filenames for easy reference.
  
**Outputs**

- RCSB_RNA_Binding_Protein.txt

- RCSB_DNA_Binding_Protein.txt

- RCSB_combined_PDB_IDs.txt

- config.json
  
**2_ecRBPome.py**

  It automates the retrieval of PDB IDs from the ecRBPome database (RNA-binding proteins in Escherichia coli). Since the database does not offer any outputs for PDB IDs, the script uses Selenium to navigate the webpage and interact with elements to extract the relevant data. It first loads the ecRBPome database page and sorts the "Cross Reference (PDB)" column twice to ensure that rows containing PDB IDs appear at the top. It then iterates through each row and extracts PDB identifiers from hyperlinks in the fifth column (Cross Reference (PDB)).
  
  The script handles pagination and also clicking the "Next" button to move through all available pages and checking for termination conditions at the same time. All PDB IDs are saved in ecRBPome_pdb_ids.txt that allow for later comparison with the RCSB dataset. This script is essential for gathering structural data on RNA-binding proteins that are missing in RCSB searched results and will later be compared with the PDB IDs retrieved from RCSB.
  
**Outputs**

- ecRBPome_pdb_ids.txt

**3_ecRBPome_dups.py**

  This file provide integrity by removing duplicate PDB IDs from the extracted ecRBPome dataset. This script processes the ecRBPome_PDB_IDs.txt file by reading each line and stripping whitespace to further store the unique PDB IDs in a set (inherently prevents duplication). The unique entries are then sorted and saved in a new output file (ecRBPome_PDB_IDs_unique.txt) to make the dataset clean for further analysis.

**Outputs**

- ecRBPome_PDB_IDs_unique.txt

**4_datasets_comp.py**

  This script compares PDB IDs retrieved from ecRBPome and RCSB and identifies common and unique entries between the two datasets. It first retrieves the de-duplicated PDB IDs from ecRBPome_PDB_IDs_unique.txt and RCSB PDB IDs from the config.json file that stores the filename of the combined RCSB dataset. It determines the intersection (common IDs), IDs unique to ecRBPome and IDs unique to RCSB by use of set operations. The results are logged into three separate files (comparison_log_common.txt, comparison_log_unique_to_ecRBPome.txt, and comparison_log_unique_to_RCSB.txt for shared PDBs, IDs found only in ecRBPome and IDs exclusive to RCSB, respectively). Moreover, a comma-separated list of all unique PDB IDs (without duplicates) is saved in all_ids_comma_delimited.txt for downstream processing.

**Outputs**

- comparison_log_common.txt
  
- comparison_log_unique_to_ecRBPome.txt
  
- comparison_log_unique_to_RCSB.txt
  
- all_ids_comma_delimited.txt

**5_PDBCIF_Download.py**
  
  It automates the downloading of PDB and CIF files from the RCSB PDB using multithreading for parallel retrieval. It first reads PDB IDs from the file all_ids_comma_delimited.txt that contains a consolidated list of unique PDB identifiers.
  
  It checks a log file (downloaded_pdbs.log) that tracks previously downloaded files to avoid redundant downloads. The script attempts to download .pdb files first and if unsuccessful it tries downloading the .cif format as an alternative. Successfully downloaded files are saved to the directory RCSB_pdb_files, and all successful downloads are logged in downloaded_pdbs.log, while failed downloads are recorded in failed_downloads.txt. Additionally, all successfully downloaded CIF files are logged separately in cif_downloads.txt for future reference. 10 parallel threads, This script significantly speeds up the retrieval process by using 10 parallel threads that makes it a good option for handling large datasets.

**Outputs**
  
- .pdb and .cif files in RCSB_pdb_files directory
  
- downloaded_pdbs.log
  
- failed_downloads.txt

- cif_downloads.txt
  
**6_DLcheck.py**
  
  It verifies the completeness and consistency of downloaded PDB files by comparing the filenames in the RCSB_pdb_files directory with the expected PDB IDs listed in all_ids_comma_delimited.txt. It first reads the expected PDB IDs from the file and then checks which PDB files are present or missing by comparing the file names (without extensions) to the expected IDs. The script identifies (1) missing PDB IDs, entries that should have been downloaded but are absent, and (2) extra files, files present in the directory but not listed in all_ids_comma_delimited.txt. If any mismatches are detected, the script prints them for manual review. This validation step is there to make sure that all necessary PDB structures have been properly retrieved before further processing and to prevent errors caused by incomplete datasets.

**7_interaction_count.py**
  
  It processes PDB and CIF files to analyze ribose-protein interactions in RNA- and DNA-binding proteins. It first loads structural data from RCSB_pdb_files and parse each file with BioPython to extract molecular chains. Using an automated classification system, the script identifies RNA, DNA and protein chains based on the first residue of each chain (hybrid chains are not considered). It then applies spatial indexing (NeighborSearch) to detect ribose atoms (C1', C2', C3', C4' and O4') that are within a 10.0 Å threshold of protein atoms to account for potential interactions. The script efficiently analyzes multiple structures simultaneously by leveraging parallel processing (ProcessPoolExecutor). The results include counts of ribose-protein interactions for each atom position and are saved in a structured CSV file (ribose_interactions.csv) with columns for PDB file name, molecule type and interaction counts per ribose position.

**Outputs**

- ribose_interactions.csv
