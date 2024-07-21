import pandas as pd
from Bio import SeqIO
from utils.utils import translate_table

# read files
alignment_file = "alignment/all_aligned.fasta"
mutations_file = "/home/hagar/UPv/viruses/refs/cmv_hotspots.csv"
UL97 = (141798, 143921)
UL54 = (78194, 81922)  # important - UL54 gene is on the complement strand!

hotspots = pd.read_csv(mutations_file)
UL97_hotspots = hotspots["UL97"].tolist()
UL54_hotspots = [int(x) for x in hotspots["UL54"].tolist() if pd.isnull(x) == False and x != 'nan']

alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
for sample, record in alignment.items():
        alignment[sample] = list(str(record.seq).upper())
reference = alignment.pop("NC_006273.2")


def get_codon_pos(aa_position, start):
    nt_pos_on_gene = aa_position + start - 1
    codon_pos = (nt_pos_on_gene-2, nt_pos_on_gene-1, nt_pos_on_gene)
    return codon_pos

def get_comp_codon_pos(aa_position, end):
    nt_pos_on_gene = end - aa_position + 1
    codon_pos = (nt_pos_on_gene+2, nt_pos_on_gene+1, nt_pos_on_gene)
    return codon_pos

def get_codon(seq, codon_pos):
    return seq[codon_pos[0]-1] + seq[codon_pos[1]-1] + seq[codon_pos[2]-1]


def complement(codon):
    comp = {"A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N",
            "-": "-"}
    return comp[codon[0]] + comp[codon[1]] + comp[codon[2]]


if __name__ == '__main__':

    for sample, seq in alignment.items():
        resist = 0  # if the sample have one of the hotspot mutations => resist = 1
        samples_mutations = []
        for hotspot in UL97_hotspots:
            codon_pos = get_codon_pos(hotspot * 3, UL97[0])
            ref_codon = get_codon(reference, codon_pos)
            seq_codon = get_codon(seq, codon_pos)
            if 'N' not in seq_codon:
                ref_aa = translate_table[ref_codon]
                seq_aa = translate_table[seq_codon] if not '-' in seq_codon else '*'
                if not seq_aa == ref_aa:
                    resist = 1
                    samples_mutations.append(sample +"\tUL97:" + ref_aa + str(hotspot) + seq_aa)
    
        for hotspot in UL54_hotspots:
            codon_pos = get_comp_codon_pos(hotspot * 3, UL54[1])
            ref_codon = complement(get_codon(reference, codon_pos))
            seq_codon = complement(get_codon(seq, codon_pos))
            if 'N' not in seq_codon:
                ref_aa = translate_table[ref_codon]
                seq_aa = translate_table[seq_codon] if not '-' in seq_codon else '*'
                if not seq_aa == ref_aa:
                    resist = 1
                    samples_mutations.append(sample +"\tUL54: " + ref_aa + str(hotspot) + seq_aa)    

        with open("CMV_resistance.txt", 'w') as f:
            if resist:
                f.write("Some samples have mutations at known resistance spots: \n" + "Sample\t AA_mut\n")
                for mut in samples_mutations:
                    f.write(mut + '\n')
            else:
                f.write("The samples don't have mutations at known resistance spot.")
    
    
