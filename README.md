Functions of this Program
Data Preparation:
-Reads a reference proteome file and stores protein sequences and IDs in a dictionary.
-Reads a mapping file (Ensembl) that links gene IDs, transcript IDs, and UniProt IDs and stores them in another dictionary.
Data Processing:
-Combines data from the proteome and mapping dictionaries to create a new file that includes Gene ID, Transcript ID, UniProt ID, and Protein length.
-Reads a housekeeping gene list and maps transcript IDs to gene IDs using the mapping file.
-Adds housekeeping gene information by checking if a gene ID is in a list of housekeeping genes, marking proteins as housekeeping (1) or not (0).
Polyx Data Integration:
-Reads polyx data (regions of repetitive amino acids) from an external file (PolyX2 scanner output) and calculates the polyx count, types, and combined length of all polyX regions in each protein.
-Combines the proteome data, housekeeping status, and polyx information into a final output file.

Input files: Proteome file, Mapping file and Housekeeping gene file (txt files) with Ensembl Gene or Transcript ID
Output: Mapped proteome, Mapped proteome with "Housekeeping yes/no", Final file with additional information about polyX regions mentioned above
