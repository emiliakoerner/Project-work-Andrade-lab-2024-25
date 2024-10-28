# This workflow processes proteomic data, integrates external annotations, and augments it with polyx region analysis for further research.
# set working directory, imports
import os
os.chdir("C:/users/wrede/desktop/masterarbeit")
from collections import defaultdict

referenceproteome = 'reference_proteome.fasta'  # reference proteome of 20.000 proteins
mappingfile = 'mapping_file.txt'        # mapping file with IDs (Transcript ID, Gene ID, UniProt ID)
mapped = 'proteome_mapped.txt'          # output file with 19183 proteins
# Code to store protein data from proteome file in a dictionary (ID and sequence)
proteomedict = defaultdict(lambda: {"sequence":""}) # Create empty dictionary to store proteome data (ID -> sequence)
with open(referenceproteome, 'r') as file:  # Read proteome file
    current_id = None       # Variables to keep track of the current UniProt ID and its sequence
    current_sequence = []
    for line in file:       # Iterate through each line of the file
        if line.startswith('>'):    # If the line starts with '>', it is a header
            if current_id:          # If there is a protein being processed already, save it to the dictionary
                proteomedict[current_id]["sequence"] = ''.join(current_sequence)    # Join sequences without a divider
            current_id = line.split('|')[1] #Extract the UniProtID from the header line
            current_sequence = []           #Reset the current sequence for the next protein
        else:
            current_sequence.append(line.strip())   # If the line is not a header, store the sequence
    if current_id:              # Save the last protein sequence after the loop finishes
        proteomedict[current_id]["sequence"] = ''.join(current_sequence)

# Code to store data from mapping file in a dictionary (Gene ID -> Transcript ID, UniProt ID)
mapdict = defaultdict(lambda: {"uniprot_id": set(), "gene_id": "","transcript_id": set()}) # Create dictionary for mapping file
with open(mappingfile, 'r') as map:
    for map_line in map:
        map_columns = map_line.strip("\n").split("\t")      # Create a list of the 4 columns
        gene_id = map_columns[0]                            # Save each column in the dictionary
        transcript_id = map_columns[1]
        sp_id = map_columns[2]
        tr_id = map_columns[3]
        mapdict[gene_id]["gene_id"] = gene_id
        mapdict[gene_id]["transcript_id"].add(transcript_id)    # If a gene has multiple transcript IDs, all of them are saved in the dictionary
        mapdict[gene_id]["uniprot_id"].add(sp_id)           # Uniprot ID is a set of SwissProt and Trembl ID
        mapdict[gene_id]["uniprot_id"].add(tr_id)

# Now we have a dictionary representing the mapping file in python
# Create a new file with UniProt ID and Sequence (length) from proteome and Gene ID + T ID from mapping file
# Proteins with no entry in the mapping file are lost here
"""with open(mapped,'w') as mapped:
    mapped.write(f"Gene ID\tTranscript ID\tUniProtID\tProtein length\n")    # Write header
    for uniprot_id, proteindata in proteomedict.items():    # Iterate through each dictionary entry aka protein in the proteome
        sequence = proteindata["sequence"]                  # Variables for Sequence and Length
        length = len(sequence)
        for gene_id, data in mapdict.items():       # Then iterate through mapping dictionary
            if uniprot_id in data["uniprot_id"]:    # If uniprot from proteome is found in mapping file, IDs are saved as variables
                transcript_id = list(data["transcript_id"])[0] # Only one transcript ID is printed! Can be changed if desired
                gene_id = data["gene_id"]
                mapped.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\n")   # IDs and length are copied in output file
                break           # How to write one line per gene_id and not per uniprot id? swap for loops?"""

with open(mapped, 'w') as mapped:
    mapped.write(f"Gene ID\tTranscript ID\tUniProtID\tProtein length\n")  # Write header
    for gene_id, data in mapdict.items():       # Iterate through each gene id in the mapping dictionary
        transcript_id = list(data["transcript_id"])[0] # Only one transcript ID is printed! Can be changed if desired
        gene_id = data["gene_id"]
        for uniprot_id, proteindata in proteomedict.items(): # Then iterate through proteome dictionary
            sequence = proteindata["sequence"]  # Variables for sequence and length
            length = len(sequence)
            if uniprot_id in data["uniprot_id"]:    # If uniprot id is in mapping file,
                mapped.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\n")   # IDs and length are copied in output file
                break

#Code for adding Gene ID to housekeeping gene list using mapping dictionary (Transcript ID -> Gene ID)
hk_input = 'Housekeeping_GenesHuman.txt'# List without Gene ID, 2833 genes
hk_output = 'Housekeeping_list.txt'     # List with Gene ID, 2790 genes (Some T IDs do not exist in mapping file)
with open(hk_input,'r') as input, open(hk_output, 'w') as output:
    output.write("Gene ID\tTranscript ID\tGene Name\tRefSeq\tCCDS ID\n") # Write header
    for line in input:                      # Go through each line of the hk gene list
        parts = line.split("\t")
        transcript_id = parts[0]            # Save each column with variables
        gene_name = parts[1]
        refseq = parts[2]
        ccds_id = parts[3]
        for gene_id, data in mapdict.items():            # Go through each line in the mapping dictionary to search for matching transcript iD
            if transcript_id in data["transcript_id"]:   # If the transcript ID is found, save gene ID as variable
                gene_id = data["gene_id"]
                output.write(f"{gene_id}\t{transcript_id}\t{gene_name}\t{refseq}\t{ccds_id}")    # Copy all variables into output file
                break           # Ends inner for-loop after match to save time

# Next step: Add a column with Housekeeping yes/no -> Read HK list, if ID is in the list, write 1, if not, write
mapped = 'proteome_mapped.txt'
mapped_hk = 'proteome_mapped_hk.txt'
with open(hk_output,'r') as hk_list:
    hk_genes = set()        # Empty set for genes listed as housekeeping
    for line in hk_list:
         hk_genes.add(line.strip().split("\t")[0])  # Add each gene ID to the set
print("proteins in housekeeping set:", len(hk_genes))
########
hk_gene_count = {}
with open(hk_output, 'r') as hk_list:
    for line in hk_list:
        gene_id = line.strip().split("\t")[0]
        if gene_id in hk_gene_count:
            hk_gene_count[gene_id] += 1  # Count duplicates
        else:
            hk_gene_count[gene_id] = 1

# Print the number of unique genes and duplicates
unique_genes = len(hk_gene_count)
duplicates = [gene for gene, count in hk_gene_count.items() if count > 1]
print(f"Total unique housekeeping genes: {unique_genes}")
print(f"Number of duplicate genes: {len(duplicates)}")
#########
hk_count = 0
with open(mapped,'r') as file, open(mapped_hk,'w') as output:
    output.write(f"Gene ID\tTranscript ID\tUniProtID\tProtein length\tHousekeeping?\n") # Write header
    next(file)  #Skip header when reading the input file
    for line in file:
        parts = line.strip().split("\t")
        gene_id = parts[0]                              # Save each column as a variable
        transcript_id = parts[1]
        uniprot_id = parts[2]
        length = parts[3]
        if gene_id in hk_genes:         # If gene is in the housekeeping set, write 1
            output.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\t1\n")
            hk_count += 1
        else:                           # If not, write 0
            output.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\t0\n")

print("Number of housekeeping genes in final file", hk_count)

# Add polyx information to the file
inputfile = 'output_polyx.txt'      # output generated by polyx scanner
outputfile = 'proteome_mapped_hk_polyx.tsv' # final output file
#create dictionary to store information for each protein ID
polyxdata = defaultdict(lambda: {"polyx_count": 0, "polyx_types": set(), "total_length": 0})
with open(inputfile) as file:
    for line in file:
        parts = line.strip().split("\t")
        protein_ID = line.split('|')[1]     # variables for each column
        polyx_type = parts[3]
        start = int(parts[1])
        end = int(parts[2])
        length = end - start + 1            # Calculating the total length of polyx regions of each region
        # Update the data for polyx count, types and length:
        polyxdata[protein_ID]["polyx_count"] += 1          # Count goes up 1 for every line with that protein ID
        polyxdata[protein_ID]["polyx_types"].add(polyx_type)    # Make a list of each polyx type
        polyxdata[protein_ID]["total_length"] += length         # Length of each polyX is added to the total length


# Code to create final document by adding polyx info to the last output file
with open(mapped_hk,'r') as input, open(outputfile,'w') as output:
    output.write(f"GeneID\tTranscriptID\tUniProtID\tLength\tHk\tPolyx_count\tPolyx_types\tPolyx_length(aa)\tPption_polyx\tCount_grouped\n") # Header
    next(input)     # Skip header when reading input file
    for line in input:      # Go through each line of the file
        parts = line.strip("\n").split("\t")    # For each line, save each column in a variable
        gene_id = parts[0]
        transcript_id = parts[1]
        uniprot_id = parts[2]
        length = parts[3]
        housekeeping = parts[4]
        if uniprot_id in polyxdata:     # If the protein has a polyx region, write polyx data into output file
            data = polyxdata[uniprot_id]
            polyx_count = data["polyx_count"]
            polyx_types = data["polyx_types"]
            polyx_types_combined = "/".join(data["polyx_types"])    # Join polyxtypes (stored as set) with / as divider
            total_length = data["total_length"]
            aa_percent = round(int(total_length)/int(length), 4)
            if polyx_count < 2:
                output.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\t{housekeeping}\t{polyx_count}\t{polyx_types_combined}\t{total_length}\t{aa_percent}\t{polyx_count}\n")
            else:
                output.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\t{housekeeping}\t{polyx_count}\t{polyx_types_combined}\t{total_length}\t{aa_percent}\t>1\n")
        else:   # If protein has no polyx region, put 0 or - instead for poly x count, type and length
            output.write(f"{gene_id}\t{transcript_id}\t{uniprot_id}\t{length}\t{housekeeping}\t0\t-\t0\t0\t0\n")
