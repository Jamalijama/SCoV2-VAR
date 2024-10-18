import pandas as pd

from FastAligner import *
from DBConstuction_1 import *
from ReadWriteFile import *

dir = './data/'
ref_file = 'SARS-CoV-2_Wuhan-Hu-1_Reference_NC_045512.2.fasta'
SNP_file = 'SNP_dict.bz2'
del_file = 'delete_dict.bz2'
in_file = 'insert_dict.bz2'
id_file = 'seq_id.csv'
all_file = 'metadata_label_noscore.csv'

_, ref = ReadRefFasta(dir + ref_file)
seq_id_lst = ReadCsvAll(dir + id_file)
# print(len(seq_id_lst))
SNP_dict = ReadBZ2(dir + SNP_file)
delete_dict = ReadBZ2(dir + del_file)
insert_dict = ReadBZ2(dir + in_file)

freq_ham_score_thresh = 0.05

df = pd.DataFrame()


def AddNew_GISAID(add_seq_lst, add_seq_id_lst):
    ref_seq = ref.upper()
    seq_id_lst_old = seq_id_lst
    # print(len(seq_id_lst))

    k = 0
    while k < len(add_seq_id_lst):
        for i in range(len(seq_id_lst_old)):
            if add_seq_id_lst[k] in seq_id_lst_old[i]:
                print('The sequence', add_seq_id_lst[k], 'is already in the database !')
                add_seq_lst.remove(add_seq_lst[k])
                add_seq_id_lst.remove(add_seq_id_lst[k])
                k -= 1
                break
        k += 1

    add_SNP_dict, add_insert_dict, add_delete_dict = FindDiff_GISAID(ref_seq, add_seq_id_lst,
                                                                     add_seq_lst)

    WriteBZ2Add(dir + SNP_file, add_SNP_dict)
    WriteBZ2Add(dir + in_file, add_insert_dict)
    WriteBZ2Add(dir + del_file, add_delete_dict)
    WriteCsvAdd(dir + id_file, add_seq_id_lst)
    print('Add ', len(add_seq_id_lst), ' sequences to the database !')


def AddNew_Align(add_seq_lst, add_seq_id_lst):
    ref_seq = ref.upper()
    reflen = len(ref_seq)
    dict_locs_dnts = ref_dnt_parse(ref_seq)
    dict_locs_pnts = ref_pnt_parse(ref_seq)

    seq_id_lst_old = seq_id_lst

    k = 0
    while k < len(add_seq_id_lst):
        for i in range(len(seq_id_lst_old)):
            if add_seq_id_lst[k] in seq_id_lst_old[i]:
                print('The sequence', add_seq_id_lst[k], 'is already in the database !')
                add_seq_lst.remove(add_seq_lst[k])
                add_seq_id_lst.remove(add_seq_id_lst[k])
                k -= 1
                break
        k += 1

    seq_id_lst_new = []
    seq_lst_new = []
    inloc_lst = []
    inNts_lst = []
    for i in range(len(add_seq_lst)):
        add_seq_lst[i] = add_seq_lst[i].upper()
        aligning = FastAlign(add_seq_lst[i], ref_seq, 'dnt', 'pnt',
                             dict_locs_dnts, dict_locs_pnts, 40)
        if len(aligning[0]) != len(ref_seq):
            print('The sequence', add_seq_id_lst[i], 'mismatch with reference sequence !')
        else:
            ham_score = ham_distance0(ref_seq, aligning[0])
            freq_ham_score = ham_score / reflen
            if freq_ham_score <= freq_ham_score_thresh:
                # print(add_seq_id_lst[i])
                seq_id_lst_new.append(add_seq_id_lst[i])
                seq_lst_new.append(aligning[0])
                inloc_lst.append(aligning[3])
                inNts_lst.append(aligning[4])
            else:
                print('The sequence', add_seq_id_lst[i], 'mismatch with reference sequence !')

    add_SNP_dict, add_insert_dict, add_delete_dict = FindDiff_Align(ref_seq, seq_id_lst_new, seq_lst_new, inloc_lst,
                                                                    inNts_lst)

    WriteBZ2Add(dir + SNP_file, add_SNP_dict)
    WriteBZ2Add(dir + in_file, add_insert_dict)
    WriteBZ2Add(dir + del_file, add_delete_dict)
    WriteCsvAdd(dir + id_file, seq_id_lst_new)
    print('Add ', len(seq_id_lst_new), ' sequences to the database !')


def IdSearch(search_seq_ids):
    ref_seq = ref.upper()

    # print('The total length of the current database is : ', len(seq_id_lst))

    search_seq_lst = []
    diff_seq_lst = []
    search_seq_id_lst = []
    for seq_id in search_seq_ids:
        flag = 0
        for i in range(len(seq_id_lst)):
            if seq_id in seq_id_lst[i]:
                seq_id = seq_id_lst[i]
                flag = 1
                break

        if flag == 0:
            print('There is no such sequence', seq_id, 'in this database!')
        else:
            search_seq_id_lst.append(seq_id)
            print(seq_id, ' is in this database!')
            # seq_lst = [char for char in ref_seq]
            final_seq = ref_seq
            diff_seq = '.' * len(ref_seq)

            if insert_dict.get(seq_id) != None:
                insert_lst = insert_dict[seq_id]
                for i in range(0, len(insert_lst), 2):
                    start = insert_lst[i]
                    seq = insert_lst[i + 1]
                    if start < 0:
                        final_seq = seq + final_seq
                        diff_seq = seq + diff_seq
                    else:
                        final_seq = final_seq[:start] + seq + final_seq[start:]
                        diff_seq = diff_seq[:start] + seq + diff_seq[start:]

            if SNP_dict.get(seq_id) != None:
                SNP_lst = SNP_dict[seq_id]
                for i in range(0, len(SNP_lst), 2):
                    start = SNP_lst[i]
                    seq = SNP_lst[i + 1]
                    final_seq = final_seq[:start] + seq + final_seq[start + len(seq):]
                    diff_seq = diff_seq[:start] + seq + diff_seq[start + len(seq):]

            if delete_dict.get(seq_id) != None:
                delete_lst = delete_dict[seq_id]
                for i in range(0, len(delete_lst), 2):
                    start = delete_lst[i]
                    length = delete_lst[i + 1]
                    final_seq = final_seq[:start] + '-' * length + final_seq[start + length:]
                    diff_seq = diff_seq[:start] + '-' * length + diff_seq[start + length:]

            string = final_seq
            # print('The sequence queried is:', string)
            search_seq_lst.append(string)
            diff_seq_lst.append(diff_seq)

    # for i in range(len(search_seq_lst)):
    #     print('the length of', search_seq_id_lst[i], 'is', len(search_seq_lst[i]))
    return search_seq_id_lst, search_seq_lst, diff_seq_lst


def GenerateSNP(ref=ref, SNP_dict=SNP_dict, seq_id_lst=seq_id_lst):
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


def SNPSearch(SNPs):
    snp_id_lst = []

    SNP_locs = [int(snp[1:-1]) for snp in SNPs]
    new_bases = [snp[-1] for snp in SNPs]
    SNPs_dict = {SNP_locs[i]: new_bases[i] for i in range(len(SNP_locs))}
    SNPs_dict = dict(sorted(SNPs_dict.items(), key=lambda x: x[0]))
    SNP_locs = list(SNPs_dict.keys())
    new_bases = list(SNPs_dict.values())

    seq_id_lst = SNP_dict.keys()
    for id in seq_id_lst:
        flags = [0 for i in range(len(SNP_locs))]
        SNP_lst = SNP_dict[id]
        for i in range(0, len(SNP_lst), 2):
            start = SNP_lst[i]
            seq = SNP_lst[i + 1]
            seqlen = len(seq)
            for j in range(len(SNP_locs)):
                if start > SNP_locs[j] or start + seqlen <= SNP_locs[j]:
                    continue
                else:
                    # print(start, start + seqlen, SNP_locs[j])
                    k = SNP_locs[j] - start
                    if seq[k] == new_bases[j]:
                        flags[j] = 1
        # print(flags)
        if 0 not in flags:
            snp_id_lst.append(id)

    search_seq_id_lst, search_seq_lst, diff_seq_lst = IdSearch(snp_id_lst)
    return search_seq_id_lst, search_seq_lst


# not use and maybe not correct
def SimilarSearch(search_seq):
    seq_len = len(search_seq)
    ref_len = len(ref)
    dis_lst = []
    for i in range(0, ref_len - seq_len + 1):
        little_ref = ref[i:i + seq_len]
        dis_lst.append(ham_distance0(search_seq, little_ref))
    min_dis = min(dis_lst)
    min_index = dis_lst.index(min_dis)
    sstart = min_index
    end = sstart + seq_len
    # print(sstart, end)

    start_lst = []
    seq_lst = []
    dis_lst = []
    for seq_id in seq_id_lst:
        final_seq = ref
        start_tmp = sstart

        if insert_dict.get(seq_id) != None:
            insert_lst = insert_dict[seq_id]
            for i in range(0, len(insert_lst), 2):
                start = insert_lst[i]
                seq = insert_lst[i + 1]
                if start < 0:
                    final_seq = seq + final_seq
                else:
                    start_tmp -= len(seq)
                    final_seq = final_seq[:start] + seq + final_seq[start:]

        if SNP_dict.get(seq_id) != None:
            SNP_lst = SNP_dict[seq_id]
            for i in range(0, len(SNP_lst), 2):
                start = SNP_lst[i]
                seq = SNP_lst[i + 1]
                final_seq = final_seq[:start] + seq + final_seq[start + len(seq):]

        if delete_dict.get(seq_id) != None:
            delete_lst = delete_dict[seq_id]
            for i in range(0, len(delete_lst), 2):
                start = delete_lst[i]
                length = delete_lst[i + 1]
                final_seq = final_seq[:start] + '-' * length + final_seq[start + length:]

        seq_lst.append(final_seq)
        start_lst.append(start_tmp)
        dis_tmp = ham_distance0(search_seq, final_seq[start_tmp:end])
        dis_lst.append(dis_tmp)
    df = pd.DataFrame()
    df['seq_id'] = seq_id_lst
    df['start'] = start_lst
    df['distance'] = dis_lst
    df.to_excel(dir + 'search_similar_result.xlsx', index=False)


def DetailInf(seq_id_lst, seq_lst):
    df_ = pd.DataFrame()
    for id in seq_id_lst:
        df1 = df[df['Accession ID'] == id]
        df_ = pd.concat([df_, df1])

    df2 = df_
    df2['seq'] = seq_lst
    return df2


def ContinentDateSearch(continent, year, month):
    # print(df.shape)
    df1 = df[df['Continent'] == continent + ' ']
    # print(df1.shape)
    df2 = df1[df1['Collection Year'] == int(year)]
    # print(df2.shape)
    df3 = df2[df2['Collection Month'] == int(month)]
    # print(df3.shape)
    ids = df3['Accession ID'].tolist()
    search_seq_id_lst, search_seq_lst, diff_seq_lst = IdSearch(ids)
    df3['seq'] = search_seq_lst
    return search_seq_id_lst, search_seq_lst, diff_seq_lst, df3


if __name__ == '__main__':
    while True:
        print('-' * 100)
        print('### Currently, the database supports the following operations:')
        print('### For add operations from fasta file, please enter: 1')
        print('### For query operations based on strain name, please enter: 2')
        print('### For query operations based on SNPs, please enter: 3')
        print('### For query operations based on continent and date, please enter: 4')
        print('### No action, exit, please enter: 5')
        num = input('###### Please select the action you want to take:')
        if num == '1':
            print('### Currently, the database supports data from GISAID or others:')
            print('### For data from GISAID, please enter: 1')
            print('### For data from others, please enter: 2')
            num1 = input('###### Please select the data source:')
            fafile = input('###### Please enter the relative path and the full name of the fasta file '
                           'for the new sequence to be added:')
            add_seq_lst, add_seq_id_lst = ReadFasta(fafile)
            if num1 == '1':
                AddNew_GISAID(add_seq_lst, add_seq_id_lst)
            else:
                AddNew_Align(add_seq_lst, add_seq_id_lst)
        elif num == '2':
            print(
                '### Currently, the database supports sequence query, variations from the reference sequence query, '
                'and sequence and annotation information query:')
            print('### For sequence query and generate into FASTA, please enter: 1')
            print('### For variations from the reference sequence query and generate into FASTA, please enter: 2')
            print('### For sequence and annotation information query and generate into XLSX, please enter: 3')
            num1 = input('###### Please select the the appropriate query:')
            seq_id_path = input(
                '###### Please enter the relative path and the full name of the file '
                'in which each row contains one serial name of the sequence:')
            search_seq_ids = []
            if seq_id_path.endswith('.txt'):
                search_seq_ids = ReadTxt(seq_id_path)
            elif seq_id_path.endswith('.csv') or seq_id_path.endswith('.xlsx') or seq_id_path.endswith('.xls'):
                search_seq_ids = ReadFirstCol(seq_id_path)
            search_seq_id_lst, search_seq_lst, diff_seq_lst = IdSearch(search_seq_ids)
            if seq_id_path.endswith('.xlsx'):
                if num1 == '1':
                    WriteFasta(seq_id_path[:-5] + '_full_seq.fasta', search_seq_id_lst, search_seq_lst)
                elif num1 == '2':
                    WriteFasta(seq_id_path[:-5] + '_diff.fasta', search_seq_id_lst, diff_seq_lst)
                elif num1 == '3':
                    df = pd.read_csv(dir + all_file)
                    df = df[df['label_0.99'] == 1]
                    df2 = DetailInf(search_seq_id_lst, search_seq_lst)
                    df2.to_excel(seq_id_path[:-5] + '_detail.xlsx', index=False)
            else:
                if num1 == '1':
                    WriteFasta(seq_id_path[:-4] + '_full_seq.fasta', search_seq_id_lst, search_seq_lst)
                elif num1 == '2':
                    WriteFasta(seq_id_path[:-4] + '_diff.fasta', search_seq_id_lst, diff_seq_lst)
                elif num1 == '3':
                    df = pd.read_csv(dir + all_file)
                    df = df[df['label_0.99'] == 1]
                    df2 = DetailInf(search_seq_id_lst, search_seq_lst)
                    df2.to_excel(seq_id_path[:-4] + '_detail.xlsx', index=False)
        elif num == '3':
            print(
                '### Currently, the database supports sequence query and annotation information query:')
            print('### For sequence query and generate into FASTA, please enter: 1')
            print('### For sequence and annotation information query and generate into XLSX, please enter: 2')
            num1 = input('###### Please select the the appropriate query:')
            path = input('###### Please enter the relative path and the full name of the file '
                         'in which each row contains one SNP:')
            search_snps = []
            if path.endswith('.txt'):
                search_snps = ReadTxt(path)
            elif path.endswith('.csv') or path.endswith('.xlsx') or path.endswith('.xls'):
                search_snps = ReadFirstCol(path)
            search_seq_id_lst, search_seq_lst = SNPSearch(search_snps)
            if path.endswith('.xlsx'):
                if num1 == '1':
                    WriteFasta(path[:-5] + '_full_seq.fasta', search_seq_id_lst, search_seq_lst)
                elif num1 == '2':
                    df = pd.read_csv(dir + all_file)
                    df = df[df['label_0.99'] == 1]
                    df2 = DetailInf(search_seq_id_lst, search_seq_lst)
                    df2.to_excel(path[:-5] + '_detail.xlsx', index=False)
            else:
                if num1 == '1':
                    WriteFasta(path[:-4] + '_full_seq.fasta', search_seq_id_lst, search_seq_lst)
                elif num1 == '2':
                    df = pd.read_csv(dir + all_file)
                    df = df[df['label_0.99'] == 1]
                    df2 = DetailInf(search_seq_id_lst, search_seq_lst)
                    df2.to_excel(path[:-4] + '_detail.xlsx', index=False)
        elif num == '4':
            print(
                '### Currently, the database supports sequence query, variations from the reference sequence query, '
                'and sequence and annotation information query:')
            print('### For sequence query and generate into FASTA, please enter: 1')
            print('### For variations from the reference sequence query and generate into FASTA, please enter: 2')
            print('### For sequence and annotation information query and generate into XLSX, please enter: 3')
            num1 = input('###### Please select the the appropriate query:')
            continent = input('###### Please enter the continent:')
            year = input('###### Please enter the year:')
            month = input('###### Please enter the month:')
            df = pd.read_csv(dir + all_file)
            df = df[df['label_0.99'] == 1]
            search_seq_id_lst, search_seq_lst, diff_seq_lst, df2 = ContinentDateSearch(continent, year, month)
            if num1 == '1':
                WriteFasta('%s%s_%s_%s_full_seq.fasta' % (dir1, continent, year, month), search_seq_id_lst,
                           search_seq_lst)
            elif num1 == '2':
                WriteFasta('%s%s_%s_%s_diff.fasta' % (dir1, continent, year, month), search_seq_id_lst, diff_seq_lst)
            elif num1 == '3':
                df2.to_excel('%s%s_%s_%s_detail.xlsx' % (dir1, continent, year, month), index=False)
        elif num == '5':
            break
