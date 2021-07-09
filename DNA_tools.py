import sys
import random

def reverse_complement(DNA_seq, type = "DNA"):
    if type == "DNA":
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
                      'a':'t', 't':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    elif type == "RNA":
        complement = {'A':'U', 'U':'A', 'G':'C', 'C':'G',
                      'a':'u', 'u':'a', 'g':'c', 'c':'g',
                      'N':'N', 'n':'n'}
    else:
        sys.stderr.write("Did not regognize sequence type: %s" % type)
        sys.exit()
    revcomp = []
    for i in range(len(DNA_seq))[::-1]:
        try:
            revcomp.append(complement[DNA_seq[i]])
        except KeyError:
            sys.stderr.write("Did not recognize base: %s" % DNA_seq[i])
            sys.exit()
    return ''.join(revcomp)

def align_primers(primer1, primer2, direction = "revcomp", min_overlap = 10):
    hit = False
    if direction == "revcomp":
        target = "-" * (len(primer1)-min_overlap) + reverse_complement(primer2) + "-" * (len(primer1)-min_overlap)
    for i in range(len(target)):
        hitscore = 0
        for j in range(len(primer1)):
            try:
                if primer1[j].upper() == target[i+j].upper():
                    hitscore +=1
                    continue
            except IndexError:
                break
            if target[i+j] == "-":
                continue
            else:
                hitscore = -100
                break
        if hitscore > min_overlap:
            hit = True
            query = "-" * i + primer1 + "-" * (len(target) - i - len(primer1))
            print_query = []
            print_target = []
            for i in range(len(query)):
                if query[i] != "-" or target[i] != "-":
                    print_query.append(query[i])
                    print_target.append(target[i])
            print(''.join(print_query))
            print(''.join(print_target))
    if hit == False:
        print("No suitable alignment found")
    return hit


def PCR(template, primer1, primer2, min_anneal_len = 14):
    template = template.upper()
    templateR = reverse_complement(template).upper()
    prime_ann1 = primer1.upper()[0-min_anneal_len:]
    prime_ann2 = primer2.upper()[0-min_anneal_len:]
    seeds = False
    if prime_ann1 in template and prime_ann2 in templateR:
        seeds = "fwd"
    elif prime_ann1 in templateR and prime_ann2 in template:
        seeds = "rev"
    if seeds == False:
        print(prime_ann1)
        print(prime_ann2)
        sys.stderr.write("couldn't find binding sites for primers\n")
        sys.exit()
    elif seeds == "fwd":
        primer2 = reverse_complement(primer2)
    elif seeds == "rev":
        temp = primer2
        primer2 = reverse_complement(primer1)
        primer1 = temp
    index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
    index_end = template.find(primer2.upper()[:min_anneal_len])
    if index_start < index_end:
        return(primer1 + template[index_start:index_end] + primer2)
    elif index_start > index_end:
        template = template[(index_start + index_end)//2:] + template[:(index_start + index_end)//2]
        index_start = template.find(primer1.upper()[0-min_anneal_len:])+min_anneal_len
        index_end = template.find(primer2.upper()[:min_anneal_len])
        return(primer1 + template[index_start:index_end] + primer2)
    
def read_fasta(fasta_file):
    with open(fasta_file, "r") as file_in:
        read_next_line = False
        seq_list = []
        for line in file_in:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                read_next_line = True
            elif read_next_line == True:
                if line.startswith(">"):
                    sys.stderr.write("There should only be one sequence in the fasta. Exiting...")
                    sys.exit()
                seq_list.append(line.strip())
        file_in.close()
    seq = ''.join(seq_list)

    return (seq_name, seq)

def assemble_fragments(list_of_subseqs,
                       min_overlap_len = 15,
                       max_overlap_len = 30,
                       product_topo = "linear"):
    #Note: presently requires that all fragments be in the same orientation,
    #in order, and all uppercase (or at least matching cases for the annealing portions)
    #Will fail if one or more provided fragments is not part of the final product
    if product_topo not in ["linear"]:
        sys.stderr.write("right now this only works for linear assemblies. Should add circular framework\n")
        sys.exi()
    junctions = []
    for i in range(len(list_of_subseqs)):
        for ii in range(len(list_of_subseqs)):
            for j in range(min_overlap_len, max_overlap_len+1):
                #for seq1, cut slice of j bases from the back end
                #for seq2, cut slice of j bases from the front end
                #if slices match, stich sequences together
                if list_of_subseqs[i][-j:] == list_of_subseqs[ii][:j]:
                    junctions.append((i, ii, j, list_of_subseqs[i][-j:].upper()))
                    
    junctions = sorted(junctions)
    final_product = ['|']
    for j in junctions:
        #sys.stderr.write('\t'.join([str(k) for k in j]) + '\n')
        if j[1] - j[0] != 1:
            sys.stderr.write("Seq fragments appear to be out of order, please rearrange\n")
            sys.exit()
        if list_of_subseqs[j[0]][:-j[2]].startswith(final_product[-1]):
            final_product.append(list_of_subseqs[j[0]][len(final_product[-1]):-j[2]])
        else:
            final_product.append(list_of_subseqs[j[0]][:-j[2]])
        final_product.append(j[3])
        #sys.stderr.write(''.join(final_product) + '\n')
    final_product.append(list_of_subseqs[junctions[-1][1]][junctions[-1][2]:])
    #sys.stderr.write(''.join(final_product) + '\n')
    #print(final_product)
    return ''.join(final_product).strip('|')

translation_table = {	
	'ACC': "T", 'ATG': "M", 'ACA': "T", 
	'ACG': "T", 'ATC': "I", 'AAC': "N", 
	'ATA': "I", 'AGG': "R", 'CCT': "P", 
	'CTC': "L", 'AGC': "S", 'AAG': "K", 
	'AGA': "R", 'CAT': "H", 'AAT': "N", 
	'ATT': "I", 'CTG': "L", 'CTA': "L", 
	'ACT': "T", 'CAC': "H", 'AAA': "K", 
	'CCG': "P", 'AGT': "S", 'CCA': "P", 
	'CAA': "Q", 'CCC': "P", 'TAT': "Y", 
	'GGT': "G", 'TGT': "C", 'CGA': "R", 
	'CAG': "Q", 'CGC': "R", 'GAT': "D", 
	'CGG': "R", 'CTT': "L", 'TGC': "C", 
	'GGG': "G", 'TAG': "*", 'GGA': "G", 
	'TAA': "*", 'GGC': "G", 'TAC': "Y", 
	'GAG': "E", 'TCG': "S", 'TTT': "F", 
	'GAC': "D", 'CGT': "R", 'GAA': "E", 
	'TCA': "S", 'GCA': "A", 'GTA': "V", 
	'GCC': "A", 'GTC': "V", 'GCG': "A", 
	'GTG': "V", 'TTC': "F", 'GTT': "V", 
	'GCT': "A", 'TTA': "L", 'TGA': "*", 
	'TTG': "L", 'TCC': "S", 'TGG': "W", 
	'TCT': "S"}

def translate(DNA_sequence, frame = 0):
	DNA_sequence = DNA_sequence.upper()
	protein_sequence_list = []
	for i in range(frame,len(DNA_sequence),3):
		protein_sequence_list.append(translation_table[DNA_sequence[i:i+3]])
	return ''.join(protein_sequence_list)

codon_synonyms = {
	'CTT': ['CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
	'ATG': [],
	'AAG': ['AAA'],
	'AAA': ['AAG'],
	'ATC': ['ATA', 'ATT'],
	'AAC': ['AAT'],
	'ATA': ['ATC', 'ATT'],
	'AGG': ['AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
	'CCT': ['CCG', 'CCA', 'CCC'],
	'CTC': ['CTG', 'CTA', 'CTT', 'TTA', 'TTG'],
	'AGC': ['AGT', 'TCG', 'TCA', 'TCC', 'TCT'],
	'ACA': ['ACC', 'ACG', 'ACT'],
	'AGA': ['AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
	'CAT': ['CAC'],
	'AAT': ['AAC'],
	'ATT': ['ATC', 'ATA'],
	'CTG': ['CTA', 'CTC', 'CTT', 'TTA', 'TTG'],
	'CTA': ['CTG', 'CTC', 'CTT', 'TTA', 'TTG'],
	'ACT': ['ACC', 'ACA', 'ACG'],
	'CAC': ['CAT'],
	'ACG': ['ACC', 'ACA', 'ACT'],
	'CAA': ['CAG'],
	'AGT': ['AGC', 'TCG', 'TCA', 'TCC', 'TCT'],
	'CAG': ['CAA'],
	'CCG': ['CCT', 'CCA', 'CCC'],
	'CCC': ['CCT', 'CCG', 'CCA'],
	'TAT': ['TAC'],
	'GGT': ['GGG', 'GGA', 'GGC'],
	'TGT': ['TGC'],
	'CGA': ['AGG', 'AGA', 'CGC', 'CGG', 'CGT'],
	'CCA': ['CCT', 'CCG', 'CCC'],
	'TCT': ['AGC', 'AGT', 'TCG', 'TCA', 'TCC'],
	'GAT': ['GAC'],
	'CGG': ['AGG', 'AGA', 'CGA', 'CGC', 'CGT'],
	'TTT': ['TTC'],
	'TGC': ['TGT'],
	'GGG': ['GGT', 'GGA', 'GGC'],
	'TAG': ['TAA', 'TGA'],
	'GGA': ['GGT', 'GGG', 'GGC'],
	'TAA': ['TAG', 'TGA'],
	'GGC': ['GGT', 'GGG', 'GGA'],
	'TAC': ['TAT'],
	'GAG': ['GAA'],
	'TCG': ['AGC', 'AGT', 'TCA', 'TCC', 'TCT'],
	'TTA': ['CTG', 'CTA', 'CTC', 'CTT', 'TTG'],
	'GAC': ['GAT'],
	'TCC': ['AGC', 'AGT', 'TCG', 'TCA', 'TCT'],
	'GAA': ['GAG'],
	'TCA': ['AGC', 'AGT', 'TCG', 'TCC', 'TCT'],
	'GCA': ['GCC', 'GCG', 'GCT'],
	'GTA': ['GTC', 'GTG', 'GTT'],
	'GCC': ['GCA', 'GCG', 'GCT'],
	'GTC': ['GTA', 'GTG', 'GTT'],
	'GCG': ['GCA', 'GCC', 'GCT'],
	'GTG': ['GTA', 'GTC', 'GTT'],
	'TTC': ['TTT'],
	'GTT': ['GTA', 'GTC', 'GTG'],
	'GCT': ['GCA', 'GCC', 'GCG'],
	'ACC': ['ACA', 'ACG', 'ACT'],
	'TGA': ['TAG', 'TAA'],
	'TTG': ['CTG', 'CTA', 'CTC', 'CTT', 'TTA'],
	'CGT': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG'],
	'TGG': [],
	'CGC': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT']}
    
   
    
   
    
   
    
   
    
