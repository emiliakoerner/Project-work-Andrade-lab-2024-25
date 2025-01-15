import os
import subprocess
from collections import defaultdict
import re
import tempfile
os.chdir("C:/Users/wrede/OneDrive - JGU/Desktop/Masterarbeit/iupred3")
"""def Proteome_dictionary(proteome): #Code to store protein data from proteome in a dict (UniprotID, gene name, sequence)
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
Proteome_dictionary("reference_proteome.fasta")

def sanitize_filename(name):
    return re.sub(r'[<>:"/\\|?*]', '_', name)

def run_iupred(uniprot_id, gene_name, sequence, output_dir): # Define function to run iupred
    uniprot_id = sanitize_filename(uniprot_id)
    gene_name = sanitize_filename(gene_name or "unknown")
    output_file = f"{uniprot_id}_{gene_name}.iupred.txt" # file name consists of uniprot id and gene name
    output_path = os.path.join(output_dir, output_file) # construct file path (file name and directory specified in function)
    # Create a temporary file to store the sequence
    with tempfile.NamedTemporaryFile(delete=False) as temp_seq_file:
        temp_seq_file.write(f">{uniprot_id}\n{sequence}\n".encode())  # Write sequence in FASTA format
        temp_seq_path = temp_seq_file.name  # Path to the temporary file
    process = subprocess.run(
        ["python", "iupred3.py", temp_seq_path, "long"],  # run IUPred3 in "long" mode
        stdout=subprocess.PIPE,     # capture iupred3 output
        stderr=subprocess.PIPE      # capture error messages
        )
    if process.returncode != 0:
        print(f"Error running IUPred3 for {uniprot_id}. Return code: {process.returncode}")
        print(f"Error details: {process.stderr}")
        return None
    with open(output_path, 'w', encoding='utf-8') as file:
        file.write(process.stdout.decode())      # put iupred3 output into our output file
    return output_path
run_iupred("P12345", "BRCA1", "MESQVAQLKMNCVWTRSAQ", "iupred_results")

min_length = 19
def process_proteins(output_dir):
    for uniprot_id, protein_data in proteomedict.items():
        sequence = protein_data["sequence"]
        if len(sequence) < min_length:
            continue
        gene_name = protein_data["gene_name"]
        output_file = run_iupred(uniprot_id, gene_name, sequence, output_dir)
process_proteins("iupred_results") #took ~7 hours
"""

def calculate_values(resultsdir):
    global disorderdict
    disorderdict = defaultdict(lambda: {"disorder": "", "disordered_aas": ""})
    for filename in os.listdir(resultsdir):
        uniprot_id = filename.split("_")[0]
        filepath = os.path.join(resultsdir, filename)
        with open(filepath, 'r') as file:
            scores = []
            for line in file:
                if line.startswith("#") or not line.strip():
                    continue
                else:
                    score = float(line.strip().split()[2])
                    scores.append(score)
            disorderscores = sum(1 for s in scores if s > 0.5)
            disorder = "yes" if disorderscores > 0 else "no"
            disorderdict[uniprot_id]["disorder"] = disorder # save in dictionary
            disorderdict[uniprot_id]["disordered_aas"] = disorderscores
resultsdir = "iupred_results"
calculate_values(resultsdir)

def Final_doc_iupred(mapped, outputfile):    #Code to create final document by adding polyx info to the last output file
    with open(mapped, 'r') as input, open(outputfile, 'w') as output:
        output.write(f"GeneID\tGenename\tUniProtID\tLength\tHk\tDisorder\tdisordered_aas\tpption_disordered\n") # Header
        next(input)     # Skip header when reading input file
        for line in input:      # Go through each line of the file
            parts = line.strip("\n").split("\t")    # For each line, save each column in a variable
            gene_id = parts[0]
            gene_name = parts[1]
            uniprot_id = parts[2]
            length = (parts[3])
            housekeeping = parts[4]
            data = disorderdict[uniprot_id]
            disorder = data["disorder"]
            if disorder == "yes":
                disordered_aas = data["disordered_aas"]
                disordered_percent = round(int(disordered_aas) / int(length), 2)
                output.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t{disorder}\t{disordered_aas}\t{disordered_percent}\n")
            else:
                output.write(f"{gene_id}\t{gene_name}\t{uniprot_id}\t{length}\t{housekeeping}\t{disorder}\t0\t0\n")

Final_doc_iupred('proteome_mapped_hk.txt', 'proteome_mapped_hk_disorder.tsv')
