import os
import re
from collections import defaultdict
import requests
import time
import random
os.chdir("C:/Users/wrede/OneDrive - JGU/Desktop/Masterarbeit/model organisms")

# Function to map C elegans (WormBase) and Drosophila (Flybase) Housekeeping gene lists to Uniprot using the
# underlying data files of the Uniprot Mapping tool.
def Housekeeping_mapping_uniprot(hk_file, hk_mapping_file, hk_list):
    with open(hk_mapping_file, 'r') as mapping:
        mapping_lines = mapping.readlines()     # Mapping file is read into python storage
    with open(hk_file, 'r') as hk_file, open(hk_list, 'w') as output:
        for line in hk_file:
            wb_id = line.strip()    # Hk genes are stored and iterated through
            for line2 in mapping_lines:
                columns = line2.strip().split("\t")
                wb = columns[2]     # WB ID or Flybase ID are always in the third column of the mapping file
                if wb_id == wb:     # If a matching entry in the mapping file is found,
                    uniprot_id = columns[0]     # The Uniprot ID is stored and written into an output file with both IDs
                    output.write(f"{wb_id}\t{uniprot_id}\n")
                    break

hk_file = "fruitfly_hk_unmapped.txt"
hk_mapping_dir = "mapping files"
hk_mapping_file = os.path.join(hk_mapping_dir, "DROME_7227_idmapping.txt")
housekeeping_dir = "housekeeping_lists"
hk_list = os.path.join(housekeeping_dir, "fruitfly_hk.txt")

#Housekeeping_mapping_uniprot(hk_file, hk_mapping_file, hk_list)
# We now have a hk gene list with both IDs.

def Process_proteomes(proteome_dir, hk_dir):
    for proteome in os.listdir(proteome_dir):       # every proteome in the directory is processed
        if proteome.endswith(".fasta"):             # if it is a fasta file
            species = proteome.replace(".fasta", "")        # the file name is saved as species
            proteome_path = os.path.join(proteome_dir, proteome)
            housekeeping_file = proteome.replace(".fasta", "_hk.txt")   # the hk gene list has the same file name but with _hk
            housekeeping_file = os.path.join(hk_dir, housekeeping_file)

            mapped = f"mapped_hk/{species}_mapped.txt"      # the output file should have the same file name but with _mapped

            proteome_dict = Proteome_dictionary(proteome_path)      # calling the function that reads proteome into a dictionary
            Map_to_hklist(proteome_dict, housekeeping_file, mapped)     # calling the function that maps proteome to hk gene list

def Proteome_dictionary(proteome): #Code to store protein data from proteome in a dict (UniprotID, gene name, sequence)
    global proteomedict
    proteomedict = defaultdict(lambda: {"gene_name": None, "sequence": ""})  #Create empty dictionary (uniprot -> gn, sequence)
    with open(proteome, 'r') as file:  # Read proteome file
        current_id = None       # Variables to keep track of the current UP ID and its sequence
        current_sequence = []
        for line in file:       # Iterate through each line of the file
            if line.startswith('>'):    # If the line starts with '>', it is a header
                if current_id:          # If there is a protein being processed already, save it to the dictionary
                    proteomedict[current_id]["sequence"] = ''.join(current_sequence)    # Join sequences without a divider
                gn_match = re.search(r"GN=(\S+)", line) # look for gene name
                uniprot_match = re.search(r">.+?\|([^\|]+)\|", line) # look for UP ID

                current_id = uniprot_match.group(1) if uniprot_match else None      # put UP ID as current ID
                gene_name = gn_match.group(1) if gn_match else None                 # store gene name too
                current_sequence = []           # Reset the current sequence for the next protein
                if current_id:
                    proteomedict[current_id]["gene_name"] = gene_name       # save gene name in the dictionary
            else:
                current_sequence.append(line.strip())   # If the line is not a header, store the sequence
        if current_id:
            proteomedict[current_id]["sequence"] = ''.join(current_sequence)  # Save the last protein sequence after loop finishes
    print(len(proteomedict), "proteins in the dictionary")
    return proteomedict

# fetching gene name synonyms from uniprot. there is only one gene name per protein in the proteomes. many hk lists have
# older gene names that are now listed as synonyms in uniprot. this function fetches the json file for each uniprot ID
# from uniprot to improve mapping.
def fetch_synonyms(uniprot_id, max_retries=5):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"   # calling the rest.uniprot page using the uniprot ID
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()     # raises potential errors
            data = response.json()      # read all data from the json file into python

            synonyms = set()

            # Extract gene names from 'genes' field
            if "genes" in data:
                for gene in data["genes"]:
                    if "orderedLocusNames" in gene:
                        for locus in gene["orderedLocusNames"]:
                            synonyms.add(locus["value"])
            # Extract additional synonyms from cross-references
            if "uniProtKBCrossReferences" in data:
                for ref in data["uniProtKBCrossReferences"]:
                    if ref["database"] in ["EnsemblPlants", "TAIR", "KEGG"]:
                        for prop in ref["properties"]:
                            synonyms.add(prop["value"])

            print(f"Fetched synonyms for {uniprot_id}: {synonyms}")
            return synonyms

        except (requests.exceptions.RequestException, ValueError) as e: # in case of an error, retry in 1 second
            retries += 1
            print(f"Failed to fetch synonyms for {uniprot_id}: {e}. Retrying {retries}/{max_retries}...")
            time.sleep(1)

    print(f"Failed to fetch synonyms for {uniprot_id} after {max_retries} retries.")
    return set()


def Map_to_hklist(proteomedict, housekeeping_list, mapped):
    global hk_genes
    hk_genes = set()  # Load housekeeping genes from file into a set
    with open(housekeeping_list, 'r') as file:
        for line in file:
            words = line.strip().split("\t")
            hk_genes.update(word.lower() for word in words)     # in lower case! case sensitivity
    with open(mapped, 'w') as mapped:
        mapped.write(f"gene_name\tuniprot_id\tlength\tHk\n")  # Write header
        hk_count = 0            # Count the hk proteins from the original file
        for uniprot_id, proteindata in proteomedict.items():  # Iterate through each dictionary entry aka protein
            sequence = proteindata["sequence"]  # Variables for Sequence, Length and gene name
            length = len(sequence)
            gene_name = proteindata["gene_name"]
            hk_status = "0"
            if gene_name:       # if protein has a gene name, check for housekeeping status
                if gene_name.lower() in hk_genes:
                    hk_status = "1" # set hk to "yes"
                    hk_count += 1
                else:
                    synonyms = fetch_synonyms(uniprot_id)   # fetch synonyms only for proteins that could not be mapped
                    if any(syn.lower() in hk_genes for syn in synonyms):    # see if the gene name is in the synonyms set for this protein
                        hk_status = "1"
                        hk_count += 1
            else:
                gene_name = "-"
            """if uniprot_id.lower() in hk_genes:   # this is for mapping via uniprot id (c elegans and drosophila)
                hk_status = "1"
                hk_count += 1 """
            mapped.write(f"{gene_name}\t{uniprot_id}\t{length}\t{hk_status}\n") # write output file
        print("Hk mapping complete:", hk_count, "hk proteins found")

#running the program
proteome_dir = "proteomes"
housekeeping_dir = "housekeeping_lists"

#Process_proteomes(proteome_dir, housekeeping_dir)

# Add polyx information to the file
def Polyx_dictionary(input):   #create dictionary to store information for each protein ID
    global polyxdata
    polyxdata = defaultdict(lambda: {"polyx_count": 0, "polyx_types": [], "polyx_lengths": [], "total_length": 0})
    with open(input) as file:
        next(file)
        for line in file:
            parts = line.strip().split("\t")
            protein_info = parts[0]  # First column, e.g., sp|Q8H102|BH128_ARATH
            protein_ID_parts = protein_info.split('|')  # Split to extract protein ID

            if len(protein_ID_parts) < 2:  # Check if the split resulted in the expected number of parts
                print(f"Skipping malformed protein ID: {protein_info}")
                continue

            protein_ID = protein_ID_parts[1]  # The second part is the protein ID (e.g., Q8H102)
            # variables for each column
            polyx_type = parts[3]
            start = int(parts[1])
            end = int(parts[2])
            length = end - start + 1            # Calculating the total length of polyx regions of each region
            # Update the data for polyx count, types and length:
            polyxdata[protein_ID]["polyx_count"] += 1          # Count goes up 1 for every line with that protein ID
            polyxdata[protein_ID]["polyx_types"].append(polyx_type)    # Make a list of each polyx type
            polyxdata[protein_ID]["polyx_lengths"].append(length)
            polyxdata[protein_ID]["total_length"] += length         # Length of each polyX is added to the total length

def Create_final_doc(mapped, outputfile):    #Code to create final document by adding polyx info to the last output file
    with open(mapped, 'r') as input, open(outputfile, 'w') as output:
        output.write(f"Genename\tUniProtID\tLength\tHk\tPolyx_count\tPolyx_types\tPolyx_lengths"
                     f"\tTotal_length\tPption_polyx\tCount_grouped\n") # Header
        next(input)     # Skip header when reading input file
        for line in input:      # Go through each line of the file
            parts = line.strip("\n").split("\t")    # For each line, save each column in a variable
            gene_name = parts[0]
            uniprot_id = parts[1]
            length = parts[2]
            housekeeping = parts[3]
            #score = (annotation_scores.get(uniprot_id, "No score"))
            if uniprot_id in polyxdata:     # If the protein has a polyx region, write polyx data into output file
                data = polyxdata[uniprot_id]
                polyx_count = data["polyx_count"]
                polyx_types_combined = "/".join(data["polyx_types"])    # Join polyxtypes (stored as set) with / as divider
                polyx_lengths = "/".join(map(str, data["polyx_lengths"]))
                total_length = data["total_length"]
                aa_percent = round(int(total_length)/int(length), 4)
                if polyx_count < 2:
                    output.write(f"{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t"
                                 f"{polyx_count}\t{polyx_types_combined}\t{polyx_lengths}\t{total_length}\t{aa_percent}\t{polyx_count}\n")
                else:
                    output.write(f"{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t"
                                 f"{polyx_count}\t{polyx_types_combined}\t{polyx_lengths}\t{total_length}\t{aa_percent}\t>1\n")
            else:   # If protein has no polyx region, put 0 or - instead for poly x count, type and length
                output.write(f"{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t0\t-\t0\t0\t0\t0\n")

# main processing function for polyx data
def Process_polyxdata(mapped_dir, polyx_dir, outputdir):
    for mappedfile in os.listdir(mapped_dir):   # process every txt file in the directory
        if mappedfile.endswith("_mapped.txt"):
            species = mappedfile.replace("_mapped.txt", "") # the species name is in the file name
            print(f"processing {species}")
            polyx_file = f"{species}_polyx.txt"
            polyx_path = os.path.join(polyx_dir, polyx_file)
            mapped_path = os.path.join(mapped_dir, mappedfile)
            output_path = f"{outputdir}/{species}_mapped_polyx.tsv"     # file paths

            Polyx_dictionary(polyx_path)        # read polyx scanner output in a dictionary
            Create_final_doc(mapped_path, output_path)      # output file
            print(f"Processed {species}: created {output_path}")

Process_polyxdata("mapped_hk", "polyx_outputs", "mapped_hk_polyx")
