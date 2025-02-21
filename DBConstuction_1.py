from ReadWriteFile import *
from FastAligner import *

dirr = './data/'
dir1 = './data/aligned/'
dir2 = './data/not_aligned/'

ref_file = 'SARS-CoV-2_Wuhan-Hu-1_Reference_NC_045512.2.fasta'
#########
# need to be modified
gisaid_file = 'xxx.fasta'
other_file = 'xxx.fasta'
SNP_file = 'SNP_dict.bz2'
del_file = 'delete_dict.bz2'
in_file = 'insert_dict.bz2'
id_file = 'seq_id.csv'

quality_ratio = 0.99


def SortAndJudge(inloc_lst, inNts_lst):
    l1, l2 = (list(t) for t in zip(*sorted(zip(inloc_lst, inNts_lst))))
    inloc_lst = l1
    inNts_lst = l2
    loc_res = []
    nt_res = []
    for i in range(len(inloc_lst)):
        if len(loc_res) == 0:
            loc_res.append([inloc_lst[i]])
            nt_res.append([inNts_lst[i]])
        elif inloc_lst[i] == inloc_lst[i - 1] + 1:
            loc_res[-1].append(inloc_lst[i])
            nt_res[-1].append(inNts_lst[i])
        else:
            loc_res.append([inloc_lst[i]])
            nt_res.append([inNts_lst[i]])
    # print(loc_res, nt_res)
    return loc_res, nt_res


def BaseFreq(seq):
    begin = 265
    end = 29674
    seq = seq[begin:end]
    count = seq.count('A') + seq.count('T') + seq.count('C') + seq.count('G')
    return count / len(seq)


def FindDiff_GISAID(ref_seq, add_seq_id_lst, add_seq_lst):
    nt_lst = ['A', 'T', 'C', 'G', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
    SNP_dict = {}
    insert_dict = {}
    delete_dict = {}

    ref_seq = ref_seq.upper()
    for i in range(len(add_seq_lst)):
        add_seq_lst[i] = add_seq_lst[i].upper()
        if BaseFreq(add_seq_lst[i]) < quality_ratio:
            # print(BaseFreq(add_seq_lst[i]))
            continue
        j = 0
        flag = [0, 0, 0]

        while j < len(ref_seq):
            # print(j)
            if add_seq_lst[i][j] != ref_seq[j]:
                tmp_str = ''

                # insertion : Stores the start site and all bases
                if ref_seq[j] == '-':
                    tmp_str += add_seq_lst[i][j]
                    if flag[0] == 0:
                        insert_dict[add_seq_id_lst[i]] = []
                        flag[0] = 1
                    k = j + 1
                    while k < len(ref_seq):
                        if add_seq_lst[i][k] != ref_seq[k] and ref_seq[k] == '-':
                            tmp_str += add_seq_lst[i][k]
                        else:
                            break
                        k += 1
                    insert_dict[add_seq_id_lst[i]].append(j)
                    insert_dict[add_seq_id_lst[i]].append(tmp_str)
                    j = k

                # deletion : Stores the start site and length
                if add_seq_lst[i][j] == '-':
                    tmp_str += '-'
                    if flag[1] == 0:
                        delete_dict[add_seq_id_lst[i]] = []
                        flag[1] = 1
                    k = j + 1
                    while k < len(ref_seq):
                        if add_seq_lst[i][k] != ref_seq[k] and add_seq_lst[i][k] == '-':
                            tmp_str += '-'
                        else:
                            break
                        k += 1
                    delete_dict[add_seq_id_lst[i]].append(j)
                    delete_dict[add_seq_id_lst[i]].append(len(tmp_str))
                    j = k

                # SNV : Stores the start site and all bases
                elif add_seq_lst[i][j] in nt_lst and ref_seq[j] in nt_lst:
                    tmp_str += add_seq_lst[i][j]
                    if flag[2] == 0:
                        SNP_dict[add_seq_id_lst[i]] = []
                        flag[2] = 1
                    k = j + 1
                    while k < len(ref_seq):
                        if add_seq_lst[i][k] != ref_seq[k] and add_seq_lst[i][k] in nt_lst and ref_seq[k] in nt_lst:
                            tmp_str += add_seq_lst[i][k]
                        else:
                            break
                        k += 1
                    SNP_dict[add_seq_id_lst[i]].append(j)
                    SNP_dict[add_seq_id_lst[i]].append(tmp_str)
                    j = k

            else:
                j += 1
            print('Add', add_seq_id_lst[i], 'to the database !')

    return SNP_dict, insert_dict, delete_dict


def FindDiff_Align(ref_seq, add_seq_id_lst, add_seq_lst, inloc_lst, inNts_lst):
    nt_lst = ['A', 'T', 'C', 'G', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
    SNP_dict = {}
    insert_dict = {}
    delete_dict = {}

    ref_seq = ref_seq.upper()
    for i in range(len(add_seq_lst)):
        add_seq_lst[i] = add_seq_lst[i].upper()
        if BaseFreq(add_seq_lst[i]) < quality_ratio:
            continue

        j = 0
        flag = [0, 0, 0]

        # insertion : Stores the start site and all bases
        if len(inloc_lst[i]) != 0:
            inloc_lst1, inNts_lst1 = SortAndJudge(inloc_lst[i], inNts_lst[i])
            if flag[0] == 0:
                insert_dict[add_seq_id_lst[i]] = []
                flag[0] = 1
            for k in range(len(inloc_lst1)):
                insert_dict[add_seq_id_lst[i]].append(inloc_lst1[k][0])
                str = ''
                for l in inNts_lst1[k]:
                    str += l
                insert_dict[add_seq_id_lst[i]].append(str)

        while j < len(ref_seq):
            # print(j)
            if add_seq_lst[i][j] != ref_seq[j]:
                tmp_str = ''

                # deletion : Stores the start site and length
                if add_seq_lst[i][j] == '-':
                    tmp_str += '-'
                    if flag[1] == 0:
                        delete_dict[add_seq_id_lst[i]] = []
                        flag[1] = 1
                    k = j + 1
                    while k < len(ref_seq):
                        if add_seq_lst[i][k] != ref_seq[k] and add_seq_lst[i][k] == '-':
                            tmp_str += '-'
                        else:
                            break
                        k += 1
                    delete_dict[add_seq_id_lst[i]].append(j)
                    delete_dict[add_seq_id_lst[i]].append(len(tmp_str))
                    j = k

                # SNV : Stores the start site and all bases
                elif add_seq_lst[i][j] in nt_lst and ref_seq[j] in nt_lst:
                    tmp_str += add_seq_lst[i][j]
                    if flag[2] == 0:
                        SNP_dict[add_seq_id_lst[i]] = []
                        flag[2] = 1
                    k = j + 1
                    while k < len(ref_seq):
                        if add_seq_lst[i][k] != ref_seq[k] and add_seq_lst[i][k] in nt_lst and ref_seq[k] in nt_lst:
                            tmp_str += add_seq_lst[i][k]
                        else:
                            break
                        k += 1
                    SNP_dict[add_seq_id_lst[i]].append(j)
                    SNP_dict[add_seq_id_lst[i]].append(tmp_str)
                    j = k

            else:
                j += 1
        print('Add', add_seq_id_lst[i], 'to the database !')

    return SNP_dict, insert_dict, delete_dict


# Process sequences from GISAID
def Comparison_GISAID(ref_file, all_file):
    ref_seq_id, ref_seq = ReadRefFasta(ref_file)
    ref_seq = ref_seq.upper()
    # print('The sequence name of the reference sequence : ', ref_seq_id)
    # print('The total length of the reference sequence is : ', len(ref_seq))

    seq_lst, seq_id_lst = ReadFasta(all_file)
    SNP_dict, insert_dict, delete_dict = FindDiff_GISAID(ref_seq, seq_id_lst, seq_lst)
    return SNP_dict, insert_dict, delete_dict, seq_id_lst


# Process sequences not from GISAID
def Comparison_Align(ref_file, all_file):
    ref_seq_id, ref_seq = ReadRefFasta(ref_file)
    ref_seq = ref_seq.upper()
    reflen = len(ref_seq)
    # print('The sequence name of the reference sequence : ', ref_seq_id)
    # print('The total length of the reference sequence is : ', len(ref_seq))
    dict_locs_dnts = ref_dnt_parse(ref_seq)
    dict_locs_pnts = ref_pnt_parse(ref_seq)

    seq_lst, seq_id_lst = ReadFasta(all_file)

    seq_id_lst_new = []
    seq_lst_new = []
    inloc_lst = []
    inNts_lst = []

    # print(ref_seq)
    print('-' * 50)
    for i in range(len(seq_lst)):
        seq_lst[i] = seq_lst[i].upper()
        freq_ham_score_lst = []
        tmp_aligning_lst = []
        for j in range(300, 600, 100):
            aligning = FastAlign(seq_lst[i], ref_seq, 'dnt', 'pnt',
                                 dict_locs_dnts, dict_locs_pnts, 40, j)
            ham_score = ham_distance0(ref_seq, aligning[0])
            freq_ham_score = ham_score / reflen
            freq_ham_score_lst.append(freq_ham_score)
            tmp_aligning_lst.append(aligning)
        min_score = min(freq_ham_score_lst)
        min_idx = freq_ham_score_lst.index(min_score)
        opt_aligning = tmp_aligning_lst[min_idx]
        mapflag = opt_aligning[0]
        if mapflag:
            seqlen = len(ref_seq)
            if len(opt_aligning[0]) != len(ref_seq):
                continue
            else:
                ham_score = ham_distance0(ref_seq, opt_aligning[0])
                freq_ham_score = ham_score / seqlen
                if freq_ham_score <= 1 - quality_ratio:
                    print(seq_id_lst[i])
                    seq_id_lst_new.append(seq_id_lst[i])
                    seq_lst_new.append(opt_aligning[0])
                    inloc_lst.append(opt_aligning[3])
                    inNts_lst.append(opt_aligning[4])

    print('-' * 50)

    SNP_dict, insert_dict, delete_dict = FindDiff_Align(ref_seq, seq_id_lst_new, seq_lst_new, inloc_lst, inNts_lst)
    return SNP_dict, insert_dict, delete_dict, seq_id_lst_new


if __name__ == '__main__':
    # Add GISAID sequences
    SNP_dict, insert_dict, delete_dict, seq_id_lst = Comparison_GISAID(dirr + ref_file, dir1 + gisaid_file)
    WriteBZ2(dirr + SNP_file, SNP_dict)
    WriteBZ2(dirr + in_file, insert_dict)
    WriteBZ2(dirr + del_file, delete_dict)
    WriteCsv(dirr + id_file, seq_id_lst)
    print('The total length of the current database is : ', len(seq_id_lst))

    # Add other sequences
    SNP_dict, insert_dict, delete_dict, seq_id_lst = Comparison_Align(dirr + ref_file, dir2 + other_file)
    WriteBZ2Add(dirr + SNP_file, SNP_dict)
    WriteBZ2Add(dirr + in_file, insert_dict)
    WriteBZ2Add(dirr + del_file, delete_dict)
    WriteCsvAdd(dirr + id_file, seq_id_lst)
    print('The total length of the current database is : ', len(seq_id_lst))
