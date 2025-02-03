import requests
import json

def search_rcsb(query, organism_filter, output_file, combined_pdb_ids):
#search RCSB API a single query and append to set
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" #query json extracted realtime by sending query from UI and tracking it in network tab of inspect
    search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": query}
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "exact_match",
                        "value": organism_filter #for E Coli
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 10000 # single query without loop
            },
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined",
            "results_content_type": ["experimental"]
        }
    }
    try:
        response = requests.post(search_url, json=search_query)
        response.raise_for_status()
        search_results = response.json()
        result_set = search_results.get("result_set", [])

        pdb_ids = [entry["identifier"] for entry in result_set] # extract and append PDBID to set
        combined_pdb_ids.update(pdb_ids)

        with open(output_file, 'w') as f: #save individual results to file
            f.write(",".join(pdb_ids))
        print(f"{len(pdb_ids)} PDB IDs for '{query}' in {output_file}")

    except Exception as e:
        print(f"Error searching for '{query}': {e}")

if __name__ == "__main__":
    search_terms = ["RNA Binding Protein", "DNA Binding Protein"]  #search terms
    organism = "Escherichia coli"
    output_files = []
    combined_pdb_ids = set()  # set to store all PDB IDs
    for term in search_terms: #search for each term
        formatted_term = term.replace(' ', '_')
        output_file = f"RCSB_{formatted_term}.txt"
        output_files.append(output_file)
        search_rcsb(term, organism, output_file, combined_pdb_ids)

    # save combined PDB IDs to a single file
    combined_output_file = "RCSB_combined_PDB_IDs.txt"
    with open(combined_output_file, 'w') as f:
        f.write(",".join(combined_pdb_ids))
    print(f"combined PDB IDs in {combined_output_file}")

    # configuration for output file names
    config_data = {
        "output_files": output_files,
        "combined_output_file": combined_output_file
    }
    with open("config.json", 'w') as config_file:
        json.dump(config_data, config_file, indent=4)
    print("configuration in config.json")
