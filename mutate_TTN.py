from align_sequences import read_single_contig_fasta, smith_waterman
ttn_gene: str = read_single_contig_fasta("TTN_RefSeq_mRNA.fasta")[1]

## TEST:
## can the aligner put a sequence in the correct part of a gene?
## can the aligner identify the operations I made?

# Define a small set of mutations in one short region of the gene. 
# Take a short region of the gene
fragment_start: int = 1000
fragment_end: int = 1100

ttn_frag: str = ttn_gene[fragment_start : fragment_end]

smith_waterman(ttn_gene, ttn_frag, 1, 1, 1, 0.5)


# Deletions (should appear as '---' in ttn_frag in alignment)
ttn_frag = ttn_frag[:95] + ttn_frag[98:] # chasrs 95, 96, 97 are deleted
ttn_frag = ttn_frag[:50] + ttn_frag[55:] # chasrs 50-54 are deleted
ttn_frag = ttn_frag[2:] # chars 0 & 1 are deleted

# Insertion (should appear as '---' in ttn_gene in alignment)
ttn_frag = ttn_frag[:55] + "cccc" + ttn_frag[55:]


# Substitutions
mutations: dict = {
    "a" : "t", 
    "c" : "a",
    "g" : "a",
    "t" : "c"
    }
ttn_frag = ttn_frag[:10] + mutations[ttn_frag[10]] + ttn_frag[11:]
ttn_frag = ttn_frag[:25] + mutations[ttn_frag[25]] + ttn_frag[26:]
ttn_frag = ttn_frag[:60] + mutations[ttn_frag[60]] + ttn_frag[61:]

smith_waterman(ttn_gene, ttn_frag, 1, 1, 1, 0.5)