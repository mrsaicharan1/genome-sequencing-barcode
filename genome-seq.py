from Bio import SeqIO, Data
import numpy as np
import matplotlib.pyplot as plt 

def codon_freq_subs(s, codon, start, end):
    freq = 0
    for index in range(len(s[start:end+1])-2):
        if s[index:index+3] == codon:
            freq += 1
    return freq
    

# print(codon_freq_subs('ACTACTACTCACTGGACG', 'ACG', 0, 17))

def get_all_codons():
    codons_list = list(Data.CodonTable.unambiguous_dna_by_id[2].forward_table.keys())
    return codons_list


def extrapolate_main():
    g_seq = ''
    for seq_record in SeqIO.parse("ex1.gbk", "genbank"):
        g_seq += str(seq_record.seq)
    
    codons = get_all_codons()
    seq_length = len(g_seq)
    barcode_dict = dict()
    for codon in codons:
        freq_list = []
        index = 0
        while(index + 999 < seq_length):
            freq = codon_freq_subs(g_seq, codon, index, index+999)
            freq_list.append(freq)
            index += 1000
        barcode_dict[codon] = freq_list

    # print(len(g_seq))
    return barcode_dict

# print(extrapolate_main())

def get_seq():
    g_seq = ''
    for seq_record in SeqIO.parse("ex1.gbk", "genbank"):
        g_seq += str(seq_record.seq)
    return g_seq

# print(codon_freq_subs(get_seq(), 'GCC', 2001, 3000))

def plot_barcodes():
    codon_dict = extrapolate_main()
    ind = 1
    codon_ind = 1
    for codon in codon_dict.keys():
        for val in codon_dict[codon]:
            X = np.array(codon_ind)
            Y = np.array(val)
    # Plotting point using sactter method.
            print(codon_ind, ind)
            color = str(val/255.)
            plt.scatter(codon_ind, ind, c=color)
            ind += 999
        ind = 1
        codon_ind += 1
    plt.gray()
    plt.show()
    # plt.scatter(x1, y1, label = "Barcodes generated")

print(plot_barcodes())