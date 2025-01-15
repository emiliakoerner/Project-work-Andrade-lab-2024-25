import os                           # set working directory, imports
import pandas as pd
import re
from collections import defaultdict
os.chdir("C:/Users/wrede/OneDrive - JGU/Desktop/Masterarbeit/ensembl_functions")
pd.set_option("display.max_columns", None)          # display all columns of data frames always

gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"] # column names
gtf_file = "Homo_sapiens.GRCh38.113.chr.gtf"    # read file
gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=gtf_columns, low_memory=False) # create dataframe
genes_df = gtf_df[gtf_df["feature"] == "gene"]                  # one for genes
transcripts_df = gtf_df[gtf_df["feature"] == "transcript"]      # one for transcripts

def parse_attributes(attributes):
    attr_dict = {}      # read "attributes" into a dictionary
    matches = re.findall(r'(\S+) "([^"]+)"', attributes)    #looks for words in quotes "" and creates a tuple for each attribute
    for key, value in matches:
        attr_dict[key] = value                              #puts each key (for example gene ID) and its value in a dictionary
    return pd.Series({                                      # creates a series with the key-value pairs
        "gene_id": attr_dict.get("gene_id", None),
        "gene_name": attr_dict.get("gene_name", None),
        "gene_biotype": attr_dict.get("gene_biotype", None),
        "transcript_id": attr_dict.get("transcript_id", None),
        "transcript_biotype": attr_dict.get("transcript_biotype", None)

    })
# Apply parsing function to attributes column for both genes and transcripts
genes_df = genes_df.join(genes_df['attributes'].apply(parse_attributes)) # Each series is added to the data frame, forming a new column
transcripts_df = transcripts_df.join(transcripts_df['attributes'].apply(parse_attributes))

# Filter for protein-coding genes and transcripts
coding_gs = genes_df[genes_df["gene_biotype"] == "protein_coding"]
coding_ts = transcripts_df[transcripts_df["transcript_biotype"] == "protein_coding"]
coding_ts = coding_ts.drop(columns=['score', 'frame', 'attributes']) # Drop unnecessary columns

# Select only the gene location information (start and end) from protein_coding_genes (other columns are not needed for now)
gene_location_info = coding_gs[['gene_id', 'start', 'end']]
merged_transcripts = coding_ts.merge(gene_location_info,            # Merge gene df into transcript df on gene ID
    on='gene_id', how='left',  # left join to keep all transcripts
    suffixes=('_transcript', '_gene'))

mapping_ensembl = "mapping_ensembl.txt"
#print(merged_transcripts.head())       # print to check
merged_transcripts.to_csv(mapping_ensembl, index=False, sep="\t")   #write df into a tsv file
# This is the finished mapping file with all the data that is relevant for now       
print("Protein-coding genes:", coding_gs.shape[0])          #check numbers of genes and transcripts
print("Protein-coding transcripts:", coding_ts.shape[0])

##################################################################################################################

def Proteome_dictionary(proteome): # Code to store protein data from proteome file in a dictionary (UniprotID, gene name and sequence)
    global proteomedict
    proteomedict = defaultdict(lambda: {"gene_name": None,"sequence": ""}) # Create empty dictionary to store proteome data (uniprot -> gn, sequence)
    with open(proteome, 'r') as file:  # Read proteome file
        current_id = None       # Variables to keep track of the current UP ID and its sequence
        current_sequence = []
        for line in file:       # Iterate through each line of the file
            if line.startswith('>'):    # If the line starts with '>', it is a header
                if current_id:          # If there is a protein being processed already, save it to the dictionary
                    proteomedict[current_id]["sequence"] = ''.join(current_sequence)    # Join sequences without a divider
                gn_match = re.search(r"GN=([^\s]+)", line) # look for gene name
                uniprot_match = re.search(r">.+?\|([^\|]+)\|", line) # look for UP ID

                current_id = uniprot_match.group(1) if uniprot_match else None      # put UP ID as current ID
                gene_name = gn_match.group(1) if gn_match else None                 # store gene name too
                current_sequence = []           # Reset the current sequence for the next protein
                if current_id:
                    proteomedict[current_id]["gene_name"] = gene_name       # save gene name in the dictionary
            else:
                current_sequence.append(line.strip())   # If the line is not a header, store the sequence
        if current_id:
            proteomedict[current_id]["sequence"] = ''.join(current_sequence)    # Save the last protein sequence after the loop finishes
    print(len(proteomedict), "proteins in the dictionary")

Proteome_dictionary('reference_proteome.fasta')

def Map_dictionary(mapping_file):
    # Code to store data from mapping file in a dictionary (Transcript ID -> Gene name, Gene ID)
    global mapdict
    mapdict = defaultdict(lambda: {"transcript_id": "", "gene_name": "", "gene_id": ""}) # Create dictionary for mapping file
    with open(mapping_file, 'r') as map:
        for map_line in map:
            map_columns = map_line.strip().split("\t")      # Create a list of the columns
            gene_id = map_columns[6]
            transcript_id = map_columns[9]
            gene_name = map_columns[7]
            mapdict[transcript_id]["gene_id"] = gene_id     # Save the value of each column in the dictionary for each line
            mapdict[transcript_id]["transcript_id"] = transcript_id
            mapdict[transcript_id]["gene_name"] = gene_name

Map_dictionary("mapping_ensembl.txt")

referenceproteome = 'reference_proteome.fasta'
housekeeping_file = 'Housekeeping_GenesHuman.txt'
housekeeping_genes = set()              # read housekeeping gene list into a set for quick access
missing = set()
with open(housekeeping_file, 'r') as hk_file, open(referenceproteome, 'r') as proteome:
    proteome_genes = set()
    for line in proteome:
        if line.startswith(">"):  # Identify header lines
            # Extract the gene name after "GN="
            parts = line.split("GN=")
            if len(parts) > 1:  # Ensure GN= is present
                gene_name = parts[1].split()[0]  # Take the first word after GN=
                proteome_genes.add(gene_name)

    next(hk_file)   #skip header
    for line in hk_file:
        genename = line.strip().split("\t")[1]  # isolate the gene names
        if genename in proteome_genes:
            housekeeping_genes.add(genename)     # isolate the gene names
        else:
            missing.add(genename)
print("unique gene names in housekeeping file:",len(housekeeping_genes)+len(missing))

###### 'missing' is read into a file, which needs to be added manually (no code yet)!

missing_hks_file = 'missing_hks.txt'   # Code for that in retired_code.py. Gene names from HRT that are not in mapping file.
# Gene names have to be added to this file before using it, manually or with code retrieving them from UniProt
missing_hks = set()         # Load missing housekeeping genes from file into a set
with open(missing_hks_file, 'r') as file:
    next(file) # skip header
    for line in file:
        missing_hks.add(line.strip().split("\t")[1])  # Strip newline and add gene to set

def Map_proteins(mapped):
    with open(mapped, 'w') as mapped:
        mapped.write(f"gene_id\tgene_name\tuniprot_id\tlength\tHk\n")  # Write header
        unmapped_proteins = 0       # count proteins that cant be mapped
        matched_genes = set()       # set to avoid duplications
        hk_genes = 0            # Count the hk proteins from the original file and
        hk_genes_from_missing_hks_file = 0      # from the additional file (manually)
        no_gn = []      # list of proteins without gene names
        unmapped = []   # list of proteins that cant be mapped
        for uniprot_id, proteindata in proteomedict.items():  # Iterate through each dictionary entry aka protein in the proteome
            sequence = proteindata["sequence"]  # Variables for Sequence, Length and gene name
            length = len(sequence)
            gene_name = proteindata["gene_name"]
            if gene_name:       # if protein has a gene name, check for housekeeping status
                if gene_name in housekeeping_genes:
                    hk_status = "1" # set hk to "yes"
                    hk_genes += 1
                elif gene_name in missing_hks:
                    hk_status = "1" # hk yes
                    hk_genes_from_missing_hks_file += 1
                else: hk_status = "0" # hk no!
                found = False       # flag for successfull mapping or not
                for transcript, mapdata in mapdict.items():       # Iterate through each gene id in the mapping dictionary
                    if gene_name in mapdata["gene_name"]:    # If gene name is in mapping dictionary, it can be mapped
                        gene_id = mapdata["gene_id"]            # retrieve gene id from mapping dict
                        mapped.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{hk_status}\n")   # IDs and length are written in output file
                        matched_genes.add(uniprot_id)           # add protein to the set to avoid duplications
                        found = True                            # mapped successfully
                        break
                if not found:
                    unmapped.append((uniprot_id, gene_name, length, hk_status))    # if not mapped, add to list of unmapped proteins
            else:
                no_gn.append((uniprot_id, "-", length))     # if there is no gene name, add it to the list for proteins without gene name
        for uniprot_id, gene_name, length in no_gn:         # add all proteins without gene name to the file
            mapped.write(f"-\t{gene_name}\t{uniprot_id}\t{length}\t{hk_status}\n")
        for uniprot_id, gene_name, length, hk_status in unmapped:      # add all unmapped proteins to the bottom of the file
            if uniprot_id not in matched_genes:
                mapped.write(f"-\t{gene_name}\t{uniprot_id}\t{length}\t{hk_status}\n")
                unmapped_proteins += 1                      # count unmapped proteins for crosschecking the numbers (add up to 20656?)
    print("from original file:", hk_genes, "from missing file:", hk_genes_from_missing_hks_file)
    print(unmapped_proteins, "proteins could not be mapped")

Map_proteins('proteome_mapped_hk.txt')


# Add polyx information to the file
def Polyx_dictionary(input, output):   #create dictionary to store information for each protein ID
    global polyxdata
    polyxdata = defaultdict(lambda: {"polyx_count": 0, "polyx_types": set(), "total_length": 0})
    with open(input) as file:
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

Polyx_dictionary('output_polyx.txt', 'proteome_mapped_hk_polyx.tsv')

def Create_final_doc(mapped, outputfile):    #Code to create final document by adding polyx info to the last output file
    with open(mapped, 'r') as input, open(outputfile, 'w') as output:
        output.write(f"GeneID\tGenename\tUniProtID\tLength\tHk\tPolyx_count\tPolyx_types\tPolyx_length(aa)\tPption_polyx\tCount_grouped\n") # Header
        next(input)     # Skip header when reading input file
        for line in input:      # Go through each line of the file
            parts = line.strip("\n").split("\t")    # For each line, save each column in a variable
            gene_id = parts[0]
            gene_name = parts[1]
            uniprot_id = parts[2]
            length = parts[3]
            housekeeping = parts[4]
            if gene_name != "-":
                if uniprot_id in polyxdata:     # If the protein has a polyx region, write polyx data into output file
                    data = polyxdata[uniprot_id]
                    polyx_count = data["polyx_count"]
                    polyx_types = data["polyx_types"]
                    polyx_types_combined = "/".join(data["polyx_types"])    # Join polyxtypes (stored as set) with / as divider
                    total_length = data["total_length"]
                    aa_percent = round(int(total_length)/int(length), 4)
                    if polyx_count < 2:
                        output.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t{polyx_count}\t{polyx_types_combined}\t{total_length}\t{aa_percent}\t{polyx_count}\n")
                    else:
                        output.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t{polyx_count}\t{polyx_types_combined}\t{total_length}\t{aa_percent}\t>1\n")
                else:   # If protein has no polyx region, put 0 or - instead for poly x count, type and length
                    output.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t0\t-\t0\t0\t0\n")

Create_final_doc('proteome_mapped_hk.txt', 'proteome_mapped_hk_polyx.tsv')
