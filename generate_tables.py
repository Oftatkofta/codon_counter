import pandas as pd
from Bio import SeqIO
import seaborn as sns
from codon_counter import *
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np



# Read data from an Excel file into a Pandas DataFrame
df = pd.read_excel("Shigella_mutants-proteomics_results-run2.xlsx")

# Load genomic sequences from two files
chromosome = next(SeqIO.parse("CP037923.1.txt", "embl"))
plasmid = next(SeqIO.parse("CP037924.1.txt", "embl"))

# Create dictionaries to store amino acid counts for 'AGA', 'UUA', 'GUA', 'GGA', 'CAA', 'AAA', and 'GAA'
r_counts = make_aa_dict(chromosome, 'AGA')
r_counts.update(make_aa_dict(plasmid, 'AGA'))

l_counts = make_aa_dict(chromosome, 'UUA')
l_counts.update(make_aa_dict(plasmid, 'UUA'))

v_counts = make_aa_dict(chromosome, 'GUA')
v_counts.update(make_aa_dict(plasmid, 'GUA'))

g_counts = make_aa_dict(chromosome, 'GGA')
g_counts.update(make_aa_dict(plasmid, 'GGA'))

q_counts = make_aa_dict(chromosome, 'CAA')
q_counts.update(make_aa_dict(plasmid, 'CAA'))

k_counts = make_aa_dict(chromosome, 'AAA')
k_counts.update(make_aa_dict(plasmid, 'AAA'))

e_counts = make_aa_dict(chromosome, 'GAA')
e_counts.update(make_aa_dict(plasmid, 'GAA'))

# Create dictionaries to store codon counts for 'AGA', 'UUA', 'GUA', 'GGA', 'CAA', 'AAA', and 'GAA'
aga_counts = make_codon_dict(chromosome, 'AGA')
aga_counts.update(make_codon_dict(plasmid, 'AGA'))

uua_counts = make_codon_dict(chromosome, 'UUA')
uua_counts.update(make_codon_dict(plasmid, 'UUA'))

gua_counts = make_codon_dict(chromosome, 'GUA')
gua_counts.update(make_codon_dict(plasmid, 'GUA'))

gga_counts = make_codon_dict(chromosome, 'GGA')
gga_counts.update(make_codon_dict(plasmid, 'GGA'))

caa_counts = make_codon_dict(chromosome, 'CAA')
caa_counts.update(make_codon_dict(plasmid, 'CAA'))

aaa_counts = make_codon_dict(chromosome, 'AAA')
aaa_counts.update(make_codon_dict(plasmid, 'AAA'))

gaa_counts = make_codon_dict(chromosome, 'GAA')
gaa_counts.update(make_codon_dict(plasmid, 'GAA'))

# Create a source dictionary to associate sequences with their source ('chromosome' or 'plasmid')
source_dict = make_source_dict(chromosome, 'chromosome')
source_dict.update(make_source_dict(plasmid, 'plasmid'))