from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict 



def codon_counter(record):
    #counts occurances of triplets in record non-pseudodene CDS
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
 

def synonymus_codon_counts(codon_count_dict, translation_table=11):
    synonymus = defaultdict(dict)
    table = CodonTable.unambiguous_rna_by_id[translation_table] 
    
    for codon, count in sorted(codon_count_dict.items()):
        aa = table[codon]
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


def count_aa_occurences(aa, peptide_sequence):
    
    return peptide_sequence.count(aa)
    

def make_protein_codon_dict(record, codon, translation_table=11):
    #returns a dict[protein_id]:(total_count_syn_codons, count_specic_codon, fraction )
    out = {}
    table = CodonTable.unambiguous_rna_by_id[translation_table] 
    aa = table.forward_table[codon]
    
    for feature in record.features:
        
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            codon_count = 0
            prot_id = feature.qualifiers['protein_id'][0]
            translation = feature.qualifiers['translation'][0]
            aa_count = translation.count(aa)
            sequence = feature.extract(record.seq)
            transcription = sequence.transcribe()
        
            for cdn in range(0, len(transcription), 3):
                cdn_seq = transcription[cdn:cdn+3]
                if cdn_seq == codon: 
                    codon_count += 1
            
            if aa_count != 0:
                fraction = codon_count/aa_count
            else:
                fraction = 0
                
            out[prot_id] = (aa_count, codon_count, fraction)
    return out
 
def make_aa_dict(record, codon, translation_table=11):
    #returns a dict[protein_id]:total_count_syn_codons
    out = {}
    table = CodonTable.unambiguous_rna_by_id[translation_table] 
    aa = table.forward_table[codon]
    
    for feature in record.features:
        
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            prot_id = feature.qualifiers['protein_id'][0]
            translation = feature.qualifiers['translation'][0]
            aa_count = translation.count(aa)

            out[prot_id] = aa_count
    return out
 
def make_codon_dict(record, codon, translation_table=11):
    #returns a dict[protein_id]: count_specic_codon
    out = {}
    table = CodonTable.unambiguous_rna_by_id[translation_table] 
    aa = table.forward_table[codon]
    
    for feature in record.features:
        
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            codon_count = 0
            prot_id = feature.qualifiers['protein_id'][0]
            sequence = feature.extract(record.seq)
            transcription = sequence.transcribe()
        
            for cdn in range(0, len(transcription), 3):
                cdn_seq = transcription[cdn:cdn+3]
                if cdn_seq == codon: 
                    codon_count += 1

            out[prot_id] = codon_count
    return out
 
def make_source_dict(record, source):
    out = {}
    for feature in record.features:
        
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            prot_id = feature.qualifiers['protein_id'][0]
            out[prot_id] = source
    return out
    
def print_codon_usage(codon_count_dict)


if __name__ == "__main__":

    chromosome =  next(SeqIO.parse("CP037923.1.txt", "embl"))
    plasmid = next(SeqIO.parse("CP037924.1.txt", "embl"))
    
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


    