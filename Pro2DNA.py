#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 16:41:08 2021

@author: bryan
"""

def main(protein_fasta_in, dna_fasta_out, email_address, force_match, run_blast):
    #Set up and read sequence from file

    Entrez.email = email_address
    query_name, query_seq = read_fasta(protein_fasta_in)
    
    #re-writing fasta should strip out any extra info that might confuse BLAST
    query_fasta = ''.join([">", query_name, '\n', query_seq])
    
    #Send search request to NCBI using default parameters for tblastn
    if run_blast == "True":
        print("waiting for BLAST to run on NCBI servers...")
        result_handle = NCBIWWW.qblast("tblastn", database="nt", sequence=query_fasta)
        print("done with BLAST")
    
        #Write results to file, so you don't have to run the search again
        with open('BLAST_results.xml','w+') as save_file:
            blast_results = result_handle.read()
            save_file.write(blast_results)
    
            result_handle.close()    
    
    #Read out the first 10 hits and see if any are exact matches
    for record in NCBIXML.parse(open("BLAST_results.xml")):
        try:
            for i in range(10):
                #Read out the relevant information from the hit
                hit_acc = record.alignments[i].accession
                hit_start = record.alignments[i].hsps[0].sbjct_start
                hit_end = record.alignments[i].hsps[0].sbjct_end
                hit_frame = record.alignments[i].hsps[0].frame[1]    
                
                #Get full sequence from Entrez (usually fast)
                Entrez_handle = Entrez.efetch(db = "nucleotide", id = hit_acc, rettype = "fasta")
                record = SeqIO.read(Entrez_handle, "fasta" ) 

                #Pull out the aligned sequence and check the translation matches
                extracted_seq = record.seq[hit_start-1:hit_end]
                
                if hit_frame < 0:
                    extracted_seq = reverse_complement(extracted_seq)
                
                if i == 0:
                    top_extracted_seq = extracted_seq
                
                if translate(extracted_seq) == query_seq: #everything is good.
                    with open(dna_fasta_out,'w+') as fasta_out:
                        fasta_out.write(">%s_top_hit_DNA\n%s\n" % (query_name, extracted_seq))
                    sys.exit()
                #else, try the next hit (kind of a shot in the dark)
                    
        #If there are no alignments, something probably failed on the server's end
        except AttributeError:
            sys.stderr.write("Blast seems to have failed. Check BLAST_results.xml for info\n")
            sys.exit()
          
        #If there are alignments but nothing matched exact, return the first match to stderr
        if force_match == "False":
            sys.stderr.write("BLAST couldn't find an exact match to your protein seq\n")
            sys.stderr.write(">Query protein\n%s\n\n" % query_seq)
            sys.stderr.write(">Top hit protein\n%s\n\n" % translate(top_extracted_seq))
            sys.stderr.write("To write the DNA seq of the top hit to file, re-run the script with '--force True --blast False'")
            sys.exit()
    
    
    
if __name__ == "__main__":

    from optparse import OptionParser
    import sys
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez, SeqIO
    from DNA_tools import reverse_complement, read_fasta, translate
    parser = OptionParser()
    parser.add_option('--input',
          '-i',
          action = 'store',
          type = 'string',
          dest = 'input_fasta',
          help = "fasta with your protein seq of ")
    parser.add_option('--output',
          '-o',
          action = 'store',
          type = 'string',
          dest = 'output_fasta',
          help = "output tsv to which normalized scores are written")
    parser.add_option('--email',
          '-e',
          action = 'store',
          type = 'string',
          dest = 'email',
          help = "Required email contact for use of Entrez servers")
    parser.add_option('--force',
          '-f',
          action = 'store',
          type = 'str',
          dest = 'force_match',
          default = 'False',
          help = "if set to true, will write the DNA sequence of the top hit even its protein sequence doesn't perfectly match your query.")
    parser.add_option('--blast',
          '-b',
          action = 'store',
          type = 'str',
          dest = 'run_blast',
          default = 'True',
          help = "if set to false, skip actually running the blast and just read from the last xml file that was written.")

    (option, args) = parser.parse_args()

    main(option.input_fasta, option.output_fasta, option.email, option.force_match, option.run_blast)