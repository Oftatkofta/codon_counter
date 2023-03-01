from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict 


chromosome =  next(SeqIO.parse("CP037923.1.txt", "embl"))
plasmid = next(SeqIO.parse("CP037924.1.txt", "embl"))


def codon_counter(record):

    codon_count = {}
    #print(record.id)
    genome = record.seq
        
    for feature in record.features:
        
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
           
            translation = feature.qualifiers['translation']
            sequence = feature.extract(genome)
            if len(sequence)%3 != 0:
                print(feature.qualifiers['locus_tag'])
                print(feature.extract(genome))
            
            transcription = sequence.transcribe()
        
            for codon in range(0, len(transcription), 3):
                codon_seq = transcription[codon:codon+3]
                if codon_seq in codon_count:
                    codon_count[codon_seq] += 1
                else:
                    codon_count[codon_seq] = 1
    return codon_count


def relative_codon_count(codon_count_dict):

    tot = 0
    for codon, count in sorted(codon_count_dict.items()):
        tot += count

    codon_frequency = {}
    for codon, count in sorted(codon_count_dict.items()):
        codon_frequency[codon] = [count/tot * 1000, count]
    
    return codon_frequency
 

def synonymus_codon_counts(codon_count_dict, translation_table):
    synonymus = defaultdict(dict)
    
    for codon, count in sorted(codon_count_dict.items()):
        aa = translation_table[codon]
        cnt = synonymus[aa]
        cnt[codon] = count
        synonymus[aa] = cnt
    
    return synonymus


def normalized_codon_usage(synonomous_dict):
    
    normdict = defaultdict(dict)
    
    for aa, codon_dict in synonomous_dict.items():
        tot = 0
        for codon, count in codon_dict.items():
            tot += count
        for codon, count in codon_dict.items():
            normdict[aa][str(codon)] = round(count/tot, 2)
    
    return normdict

table_11 = CodonTable.unambiguous_rna_by_id[11] #Table 11 Bacterial, Archaeal, Plant Plastid
bact_table = table_11.forward_table
for codon in table_11.stop_codons:
    bact_table[str(codon)] = 'STOP'

chromosome_codons = codon_counter(chromosome)
plasmid_codons = codon_counter(plasmid)

rel_chromosome = relative_codon_count(chromosome_codons)
rel_plasmid = relative_codon_count(plasmid_codons)

chromosome_synonymus = synonymus_codon_counts(chromosome_codons, bact_table)
plasmid_synonymus = synonymus_codon_counts(plasmid_codons, bact_table)

norm_chromosome = normalized_codon_usage(chromosome_synonymus)
norm_plasmid = normalized_codon_usage(plasmid_synonymus)

for aa in sorted(norm_chromosome.keys()):
    print(aa, norm_chromosome[aa])
    print(aa, norm_plasmid[aa])


    