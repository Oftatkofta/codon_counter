from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict 



def codon_counter(record):
    """
    Count the occurrences of codons (triplets) in a genomic record's non-pseudogene coding sequences (CDS).

    Parameters:
    record: Bio.SeqRecord
        A genomic record containing DNA sequence and feature information.

    Returns:
    dict
        A dictionary where keys are codons (triplets) and values are the counts of each codon in the CDS.

    This function iterates through the features in the genomic record and counts codon occurrences within
    non-pseudogene CDS regions. It returns a dictionary with codon counts.
    """

    # Initialize a dictionary to store codon counts
    codon_count = {}

    # Get the genome sequence from the input record
    genome = record.seq

    # Iterate through the features in the genomic record
    for feature in record.features:
        # Check if the feature is a coding sequence (CDS) and not a pseudogene
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            # Get the translation of the CDS (if available)
            translation = feature.qualifiers.get('translation', '')

            # Extract the DNA sequence of the CDS from the genome
            sequence = feature.extract(genome)

            # Check if the length of the sequence is divisible by 3 (valid codons)
            if len(sequence) % 3 != 0:
                print("Invalid CDS:")
                print("Locus Tag:", feature.qualifiers.get('locus_tag', 'N/A'))
                print("Sequence:", feature.extract(genome))

            # Transcribe the DNA sequence to obtain the mRNA sequence
            transcription = sequence.transcribe()

            # Iterate through the mRNA sequence in steps of 3 to extract codons
            for codon in range(0, len(transcription), 3):
                codon_seq = transcription[codon:codon + 3]

                # Count the occurrence of each codon and store it in the codon_count dictionary
                if codon_seq in codon_count:
                    codon_count[codon_seq] += 1
                else:
                    codon_count[codon_seq] = 1

    # Return the codon count dictionary
    return codon_count




def relative_codon_count(codon_count_dict):
    """
    Calculate the relative frequency of each codon in a given dictionary of codon counts.

    Parameters:
    codon_count_dict (dict): A dictionary where keys are codons (triplet nucleotide sequences) 
                            and values are the corresponding counts.

    Returns:
    dict: A dictionary where keys are codons, and values are lists containing:
        - Relative frequency per 1000 codons (float)
        - Count of the codon in the input (int)
    """

    # Calculate the total count of all codons
    tot = 0
    for codon, count in sorted(codon_count_dict.items()):
        tot += count

    # Calculate the relative frequency of each codon and store it in a new dictionary
    codon_frequency = {}
    for codon, count in sorted(codon_count_dict.items()):
        relative_freq = (count / tot) * 1000
        codon_frequency[codon] = [relative_freq, count]

    return codon_frequency

 


def synonymous_codon_counts(codon_count_dict, translation_table_nr=11):
    """
    Group codon counts by amino acid, considering synonymous codons.

    Parameters:
    codon_count_dict (dict): A dictionary where keys are codons (triplet nucleotide sequences) 
                            and values are the corresponding counts.
    translation_table_nr (int): The translation table identifier used for amino acid mapping.
                            Default is 11, which corresponds to the standard genetic code.

    Returns:
    dict: A nested dictionary where:
        - Outer keys are amino acids (single-letter code).
        - Inner dictionaries have keys as synonymous codons and values as their counts.

    This function groups codon counts by amino acid, taking into account synonymous codons.
    It returns a nested dictionary representing codon counts at the amino acid level.
    """

    # Initialize a nested dictionary to store synonymous codon counts by amino acid
    synonymous = defaultdict(dict)

    # Get the CodonTable for the specified translation table
    table = CodonTable.unambiguous_rna_by_id[translation_table_nr] 

    # Iterate through the codon counts, sorting by codon
    for codon, count in sorted(codon_count_dict.items()):
        # Get the amino acid corresponding to the codon
        aa = table[codon]

        # Retrieve the inner dictionary for the amino acid or create a new one
        codon_counts_for_aa = synonymous[aa]

        # Update the inner dictionary with the codon count
        codon_counts_for_aa[codon] = count

    return synonymous




def normalized_codon_usage(synonymous_dict):
    """
    Calculate normalized codon usage for each amino acid.

    Parameters:
    synonymous_dict (dict): A nested dictionary where:
                            - Outer keys are amino acids (single-letter code).
                            - Inner dictionaries have keys as synonymous codons and values as their counts.

    Returns:
    dict: A nested dictionary where:
        - Outer keys are amino acids (single-letter code).
        - Inner dictionaries have keys as synonymous codons (as strings) and values as normalized usage frequencies.

    This function calculates the normalized codon usage for each amino acid based on synonymous codon counts.
    It returns a nested dictionary representing normalized codon usage frequencies.
    """

    # Initialize a nested dictionary to store normalized codon usage
    normdict = defaultdict(dict)

    # Iterate through amino acids and their corresponding codon counts
    for aa, codon_dict in synonymous_dict.items():
        # Calculate the total count of synonymous codons for the amino acid
        total_count = sum(codon_dict.values())

        # Iterate through synonymous codons and calculate normalized frequencies
        for codon, count in codon_dict.items():
            # Calculate the normalized frequency and round it to two decimal places
            normalized_freq = round(count / total_count, 2)

            # Store the normalized frequency in the nested dictionary
            normdict[aa][str(codon)] = normalized_freq

    return normdict



def count_aa_occurrences(aa, peptide_sequence):
    """
    Count the occurrences of a specific amino acid in a peptide sequence.

    Parameters:
    aa (str): The single-letter code of the amino acid to count.
    peptide_sequence (str): The peptide sequence in which to count amino acid occurrences.

    Returns:
    int: The count of occurrences of the specified amino acid in the peptide sequence.

    This function takes a single-letter amino acid code (e.g., 'A' for alanine) and a peptide sequence,
    and it counts how many times the specified amino acid appears in the sequence. It returns the count as an integer.
    """
    
    # Use the count method to count the occurrences of the specified amino acid
    count = peptide_sequence.count(aa)
    
    return count

    

def make_protein_codon_dict(record, codon, translation_table=11):
    """
    Create a dictionary of protein-specific codon information.

    Parameters:
    record: Bio.SeqRecord
        A genomic record containing DNA sequence and feature information.
    codon (str): The codon (triplet nucleotide sequence) to analyze.
    translation_table (int): The translation table identifier used for amino acid mapping.
                            Default is 11, which corresponds to the standard genetic code.

    Returns:
    dict: A dictionary where keys are protein IDs, and values are tuples containing:
        - Total count of synonymous codons in the protein
        - Count of the specified codon in the protein
        - Fraction of the specified codon in the protein

    This function analyzes a genomic record and extracts information about the specified codon within
    non-pseudogene coding sequences (CDS). It returns a dictionary with protein-specific codon information.
    """

    # Initialize an empty dictionary to store protein-specific codon information
    out = {}

    # Get the CodonTable for the specified translation table
    table = CodonTable.unambiguous_rna_by_id[translation_table]

    # Determine the corresponding amino acid for the specified codon
    aa = table.forward_table[codon]

    # Iterate through the features in the genomic record
    for feature in record.features:
        # Check if the feature is a coding sequence (CDS) and not a pseudogene
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            # Initialize codon count and protein ID
            codon_count = 0
            prot_id = feature.qualifiers['protein_id'][0]
            translation = feature.qualifiers['translation'][0]

            # Count the occurrences of the specified amino acid in the protein's translation
            aa_count = translation.count(aa)

            # Extract the DNA sequence of the CDS
            sequence = feature.extract(record.seq)

            # Transcribe the DNA sequence to obtain the mRNA sequence
            transcription = sequence.transcribe()

            # Iterate through the mRNA sequence in steps of 3 to analyze codons
            for cdn in range(0, len(transcription), 3):
                cdn_seq = transcription[cdn:cdn + 3]

                # Check if the current codon matches the specified codon
                if cdn_seq == codon:
                    codon_count += 1

            # Calculate the fraction of the specified codon in the protein
            if aa_count != 0:
                fraction = codon_count / aa_count
            else:
                fraction = 0

            # Store protein-specific codon information in the output dictionary
            out[prot_id] = (aa_count, codon_count, fraction)

    return out

 
def make_aa_dict(record, codon, translation_table=11):
    """
    Create a dictionary of protein-specific counts of a specified amino acid.

    Parameters:
    record: Bio.SeqRecord
        A genomic record containing DNA sequence and feature information.
    codon (str): The codon (triplet nucleotide sequence) to analyze.
    translation_table (int): The translation table identifier used for amino acid mapping.
                            Default is 11, which corresponds to the standard genetic code.

    Returns:
    dict: A dictionary where keys are protein IDs, and values are the total counts of the specified amino acid.

    This function analyzes a genomic record and extracts information about the specified amino acid within
    non-pseudogene coding sequences (CDS). It returns a dictionary with protein-specific counts of the amino acid.
    """

    # Initialize an empty dictionary to store protein-specific amino acid counts
    out = {}

    # Get the CodonTable for the specified translation table
    table = CodonTable.unambiguous_rna_by_id[translation_table]

    # Determine the corresponding amino acid for the specified codon
    aa = table.forward_table[codon]

    # Iterate through the features in the genomic record
    for feature in record.features:
        # Check if the feature is a coding sequence (CDS) and not a pseudogene
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            # Get the protein ID and translation of the CDS
            prot_id = feature.qualifiers['protein_id'][0]
            translation = feature.qualifiers['translation'][0]

            # Count the occurrences of the specified amino acid in the protein's translation
            aa_count = translation.count(aa)

            # Store the total count of the specified amino acid in the output dictionary
            out[prot_id] = aa_count

    return out

 
def make_codon_dict(record, codon, translation_table=11):
    """
    Create a dictionary of protein-specific counts of a specified codon.

    Parameters:
    record: Bio.SeqRecord
        A genomic record containing DNA sequence and feature information.
    codon (str): The codon (triplet nucleotide sequence) to analyze.
    translation_table (int): The translation table identifier used for amino acid mapping.
                            Default is 11, which corresponds to the standard genetic code.

    Returns:
    dict: A dictionary where keys are protein IDs, and values are the counts of the specified codon.

    This function analyzes a genomic record and extracts information about the specified codon within
    non-pseudogene coding sequences (CDS). It returns a dictionary with protein-specific counts of the codon.
    """

    # Initialize an empty dictionary to store protein-specific codon counts
    out = {}

    # Get the CodonTable for the specified translation table
    table = CodonTable.unambiguous_rna_by_id[translation_table]

    # Determine the corresponding amino acid for the specified codon
    aa = table.forward_table[codon]

    # Iterate through the features in the genomic record
    for feature in record.features:
        # Check if the feature is a coding sequence (CDS) and not a pseudogene
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            # Initialize codon count and protein ID
            codon_count = 0
            prot_id = feature.qualifiers['protein_id'][0]

            # Extract the DNA sequence of the CDS and transcribe it to obtain mRNA
            sequence = feature.extract(record.seq)
            transcription = sequence.transcribe()

            # Iterate through the mRNA sequence in steps of 3 to analyze codons
            for cdn in range(0, len(transcription), 3):
                cdn_seq = transcription[cdn:cdn + 3]

                # Check if the current codon matches the specified codon
                if cdn_seq == codon:
                    codon_count += 1

            # Store the count of the specified codon in the output dictionary
            out[prot_id] = codon_count

    return out

 
def make_source_dict(record, source):
    """
    Helper function to distinguish chromosomal and plasmid genes
    
    Creates a dictionary mapping protein IDs to a specified source identifier.

    Parameters:
    record: Bio.SeqRecord
        A genomic record containing DNA sequence and feature information.
    source (str): The source identifier to be associated with protein IDs.

    Returns:
    dict: A dictionary where keys are protein IDs, and values are the specified source identifier.

    This function iterates through the features in the genomic record, specifically focusing on coding sequences (CDS)
    that are not pseudogenes. It associates protein IDs with the provided source identifier and returns a dictionary.
    """
    
    # Initialize an empty dictionary to store protein IDs mapped to the source
    out = {}

    # Iterate through the features in the genomic record
    for feature in record.features:
        # Check if the feature is a coding sequence (CDS) and not a pseudogene
        if (feature.type == "CDS") and not (feature.qualifiers.get('pseudo', False)):
            # Get the protein ID associated with the CDS
            prot_id = feature.qualifiers['protein_id'][0]

            # Map the protein ID to the provided source identifier
            out[prot_id] = source

    return out



    




    