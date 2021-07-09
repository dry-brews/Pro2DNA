# Pro2DNA

Python script to search for a DNA sequence that encodes your protein sequence of interest. You provide a fasta file with a protein sequence and you get a fasta file with the DNA sequence that encodes that protein. The script uses Biopython to do a tBLASTn search of the NCBI nonredundant database, then returns the top hit if it is an exact match to your protein.

Some important notes:
- If you use the script as part of a larger shell script, do not query the NCBI databases too frequently. They will be upset if you send them more than one request every few seconds.
- If BLAST doesn't find an exact match for your protein, you can tell it to write the DNA to file anyway using '--force True'
- You must provide an email address when querying NCBI. You are unlikely to actually receive an email from them.

usage: python Pro2DNA.py -i <my_protein_seq.fasta> -o <my_DNA_seq.fasta>

The only notable dependency is biopython, which can be installed with \\
pip install Bio \\
To run Pro2DNA, you must have DNA_tools.py in the same directory
