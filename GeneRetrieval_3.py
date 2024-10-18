from ReadWriteFile import *

dir = './data/'
ref_file = 'SARS-CoV-2_Wuhan-Hu-1_Reference_NC_045512.2.fasta'
SNP_file = 'SNP_dict.bz2'
del_file = 'delete_dict.bz2'
in_file = 'insert_dict.bz2'
id_file = 'seq_id.csv'
all_file = 'metadata_label_noscore.csv'
sampled_file = '2024-02-04_unmasked_SARS-COV-2_sampled_200.csv'

_, ref = ReadRefFasta(dir + ref_file)
seq_id_lst = ReadCsvAll(dir + id_file)
# print(len(seq_id_lst))
SNP_dict = ReadBZ2(dir + SNP_file)
delete_dict = ReadBZ2(dir + del_file)
insert_dict = ReadBZ2(dir + in_file)

seg_name_lst = ['ORF1ab_1', 'ORF1ab_2', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF10']
start_loc_lst = [265, 13467, 21562, 25392, 26244, 26522, 27201, 27393, 27755, 27893, 28273, 29557]
end_loc_lst = [13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674]
seg_number_lst = [0 for _ in range(len(seg_name_lst))]
len_lst = [end_loc_lst[i] - start_loc_lst[i] for i in range(len(start_loc_lst))]

amino_list = ['F', 'F',
              'L', 'L', 'L', 'L', 'L', 'L',
              'S', 'S', 'S', 'S', 'S', 'S',
              'Y', 'Y',
              '*', '*',
              'C', 'C',
              '*',
              'W',
              'P', 'P', 'P', 'P',
              'H', 'H',
              'Q', 'Q',
              'R', 'R', 'R', 'R', 'R', 'R',
              'I', 'I', 'I',
              'M',
              'T', 'T', 'T', 'T',
              'N', 'N',
              'K', 'K',
              'V', 'V', 'V', 'V',
              'A', 'A', 'A', 'A',
              'D', 'D',
              'E', 'E',
              'G', 'G', 'G', 'G']
codon_table = ['TTT', 'TTC',
               'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
               'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',
               'TAT', 'TAC',
               'TAA', 'TAG',
               'TGT', 'TGC',
               'TGA',
               'TGG',
               'CCT', 'CCC', 'CCA', 'CCG',
               'CAT', 'CAC',
               'CAA', 'CAG',
               'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'ATT', 'ATC', 'ATA',
               'ATG',
               'ACT', 'ACC', 'ACA', 'ACG',
               'AAT', 'AAC',
               'AAA', 'AAG',
               'GTT', 'GTC', 'GTA', 'GTG',
               'GCT', 'GCC', 'GCA', 'GCG',
               'GAT', 'GAC',
               'GAA', 'GAG',
               'GGT', 'GGC', 'GGA', 'GGG']
codon_dict = dict(zip(codon_table, amino_list))


# generate SNVs for seq_id_lst
def GenerateSNV(ref=ref, SNP_dict=SNP_dict, seq_id_lst=seq_id_lst):
    snp_lst = []
    for id in seq_id_lst:
        id_snp_lst = []
        if SNP_dict.get(id) != None:
            SNP_lst = SNP_dict[id]
            for i in range(0, len(SNP_lst), 2):
                start = SNP_lst[i]
                seq = SNP_lst[i + 1]
                for j in range(len(seq)):
                    id_snp_lst.append(ref[start + j] + str(start + j) + seq[j])
        snp_lst.append(id_snp_lst)
    return snp_lst


def cds_translator(seq):
    seqlen = len(seq)
    Protein = ''
    for codon_i in range(0, seqlen, 3):
        codon = seq[codon_i:codon_i + 3]
        if codon in codon_table:
            Protein = Protein + codon_dict[codon]
        else:
            Protein = Protein + '~'
    return Protein[:-1]


# Generate one CDS for search_seq_id_lst
def GenerateCDS(cds, search_seq_id_lst, begin=21562, end=25384):
    cds_seq_lst = []
    protein_seq_lst = []
    for seq_id in search_seq_id_lst:
        begin_1 = begin
        end_1 = end
        final_S = ref
        if insert_dict.get(seq_id) != None:
            insert_lst = insert_dict[seq_id]
            for i in range(0, len(insert_lst), 2):
                start = insert_lst[i]
                seq = insert_lst[i + 1]
                if start < 0:
                    final_S = seq + final_S
                else:
                    final_S = final_S[:start] + seq + final_S[start:]
                if start <= begin:
                    begin_1 += len(seq)
                    end_1 += len(seq)
                if start > begin and start < end:
                    end_1 += len(seq)

        if SNP_dict.get(seq_id) != None:
            SNP_lst = SNP_dict[seq_id]
            for i in range(0, len(SNP_lst), 2):
                start = SNP_lst[i]
                seq = SNP_lst[i + 1]
                final_S = final_S[:start] + seq + final_S[start + len(seq):]

        if delete_dict.get(seq_id) != None:
            delete_lst = delete_dict[seq_id]
            for i in range(0, len(delete_lst), 2):
                start = delete_lst[i]
                length = delete_lst[i + 1]
                final_S = final_S[:start] + '-' * length + final_S[start + length:]

        string = final_S[begin_1: end_1]
        if len(string) % 3 != 0:
            print(string[:3], string[-3:])
            print(len(string))
        protein_seq_lst.append(cds_translator(string))
        cds_seq_lst.append(string)

    WriteFasta(dir + cds + '.fasta', seq_id_lst, cds_seq_lst)
    WriteFasta(dir + cds + '_protein.fasta', search_seq_id_lst, protein_seq_lst)


if __name__ == '__main__':
    # a sample of generate Spike sequence and protein sequence
    # researchers can modify cds and its location
    df_sampled = pd.read_csv(dir + sampled_file)
    search_seq_id_lst = df_sampled['Accession ID'].tolist()
    GenerateCDS('Spike', search_seq_id_lst, 21562, 25384)

    # a sample of generate SNVs in search_seq_id_lst
    snp_lst = GenerateSNV(ref, SNP_dict, search_seq_id_lst)
