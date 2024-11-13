# set working directory, imports
import os
import pandas as pd
import re
from collections import defaultdict
os.chdir("C:/Users/wrede/OneDrive - JGU/Desktop/Masterarbeit/parse_ensembl")
pd.set_option("display.max_columns", None)

gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"] # column names
gtf_file = "Homo_sapiens.GRCh38.113.chr.gtf"    # read file
gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=gtf_columns, low_memory=False) # create dataframe
genes_df = gtf_df[gtf_df["feature"] == "gene"]                  # one for genes
transcripts_df = gtf_df[gtf_df["feature"] == "transcript"]      # one for transcripts

def parse_attributes(attributes):
    attr_dict = {}      # read "attributes" into a dictionary
    matches = re.findall(r'(\S+) "([^"]+)"', attributes)    #########
    for key, value in matches:
        attr_dict[key] = value
        # Return attributes of interest, with default values if not present
    return pd.Series({
        "gene_id": attr_dict.get("gene_id", None),
        "gene_name": attr_dict.get("gene_name", None),
        "gene_biotype": attr_dict.get("gene_biotype", None),
        "transcript_id": attr_dict.get("transcript_id", None),
        "transcript_biotype": attr_dict.get("transcript_biotype", None)

    })
# Apply parsing function to attributes column for both genes and transcripts
genes_df = genes_df.join(genes_df['attributes'].apply(parse_attributes))
transcripts_df = transcripts_df.join(transcripts_df['attributes'].apply(parse_attributes))

# Filter for protein-coding genes and transcripts
coding_gs = genes_df[genes_df["gene_biotype"] == "protein_coding"]
coding_ts = transcripts_df[transcripts_df["transcript_biotype"] == "protein_coding"]
coding_ts = coding_ts.drop(columns=['attributes'])  # Drop the original attributes column:
coding_ts = coding_ts.drop(columns=['score'])
coding_ts = coding_ts.drop(columns=['frame'])

# Select only the gene location information (start and end) from protein_coding_genes
gene_location_info = coding_gs[['gene_id', 'start', 'end']]
# Merge transcripts with gene location information on gene_id
merged_transcripts = coding_ts.merge(gene_location_info,
    on='gene_id', how='left',  # left join to keep all transcripts
    suffixes=('_transcript', '_gene'))

mapping_ensembl = "mapping_ensembl.txt"
print(merged_transcripts.head())
merged_transcripts.to_csv(mapping_ensembl, index=False, sep="\t")

print("Total entries in gtf_df:", gtf_df.shape[0])
print("Total entries in genes_df:", genes_df.shape[0])
print("Total entries in transcript_df:", transcripts_df.shape[0])
print("Protein-coding genes:", coding_gs.shape[0])
print("Protein-coding transcripts:", coding_ts.shape[0])
