# Functions that read and write fasta, bz2, csv, txt, gb files
import bz2
import pickle
import pandas as pd
from Bio import SeqIO


# fasta
# read fasta
def ReadRefFasta(file):
    ref_seq = ''
    with open(file, "r") as f:
        for line in f.readlines():
            line = line.strip("\n")
            if line[:1] == ">":
                ref_id = line[1:]
            else:
                ref_seq += line
    return ref_id, ref_seq


def ReadFasta(all_file):
    seq_id_lst = []
    seq_lst = []
    seq = ''
    flag = 0
    with open(all_file, "r") as f:
        for line in f.readlines():
            line = line.strip("\n")
            if line[:1] == ">":
                if flag == 1:
                    seq_lst.append(seq)
                    seq = ''
                    flag = 0
                seq_id_lst.append(line[1:])
            else:
                seq += line
                flag = 1
        seq_lst.append(seq)
    return seq_lst, seq_id_lst


# write fasta
def WriteFasta(file, seq_id_lst, seq_lst):
    fp = open(file, 'w')
    for i in range(len(seq_id_lst)):
        print('>' + seq_id_lst[i] + '\n' + seq_lst[i], file=fp)
        # print(cds_lst[j])
    fp.close()


# write fasta for addition
def WriteFastaAdd(file, seq_id_lst, seq_lst):
    fp = open(file, 'a')
    for i in range(len(seq_id_lst)):
        print('>' + seq_id_lst[i] + '\n' + seq_lst[i], file=fp)
        # print(cds_lst[j])
    fp.close()


# write fasta from csv
def CsvToFasta(file1, file2):
    df = pd.read_csv(file1)
    fp = open(file2, 'w')
    for i in range(df.shape[0]):
        print('>' + df['seq_id'][i] + '\n' + df['seq_degen'][i], file=fp)
        # print(cds_lst[j])
    fp.close()


# bz2
def ReadBZ2(file):
    data = {}
    with bz2.BZ2File(file, "rb") as f:
        while True:
            try:
                my_dict = pickle.load(f)
                data = {**data, **my_dict}
            except EOFError:
                return data


def WriteBZ2(file, my_dict):
    with bz2.open(file, 'wb') as f:
        pickle.dump(my_dict, f)


def WriteBZ2Add(file, my_dict):
    with bz2.BZ2File(file, 'ab') as file:
        pickle.dump(my_dict, file)


# csv
def WriteCsv(file, lst):
    df = pd.DataFrame()
    num_lst = [num + 1 for num in range(len(lst))]
    df['num'] = num_lst
    df['seq_id'] = lst
    df.to_csv(file, index=False)


def WriteCsvAdd(file, lst):
    df = pd.read_csv(file)
    seq_id_lst = df['seq_id'].tolist()
    for item in lst:
        seq_id_lst.append(item)
    WriteCsv(file, seq_id_lst)


def ReadFirstCol(file):
    df = pd.DataFrame()
    if file.endswith('.csv'):
        df = pd.read_csv(file)
    elif file.endswith('.xlsx') or file.endswith('.xls'):
        df = pd.read_excel(file)
    lst = df.iloc[:, 0]
    return lst


def ReadCsvAll(file):
    df = pd.read_csv(file)
    seq_id_lst = df['seq_id'].tolist()
    return seq_id_lst


def ReadCsvSome(file, num_lst):
    df = pd.read_csv(file)
    seq_id_lst = df['seq_id'].tolist()
    seq_id_lst_num = []
    for num in num_lst:
        seq_id_lst_num.append(seq_id_lst[num - 1])
    return seq_id_lst_num


# txt
def ReadTxt(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        lst = [line.strip() for line in lines]
        return lst


def WriteTxt(file, seq_id_lst):
    with open(file, 'w') as f:
        for id_ in seq_id_lst:
            f.write(id_)


def WriteTxtAdd(file, seq_id_lst):
    with open(file, 'a') as f:
        for id_ in seq_id_lst:
            f.write(id_)


# Read gb and write to fasta and csv
def ReadGbAndWrite(path, file_obj):
    seq_id_lst = []
    strain_name_lst = []
    seq_lst = []
    seq_len_lst = []
    date_lst = []
    country_lst = []
    for record in SeqIO.parse(path + file_obj, 'genbank'):
        print(record.description)
        seq_id_lst.append(record.id + ' ' + record.description)
        for feature in record.features:
            feature_type = feature.type
            qualifiers = feature.qualifiers
            # print(feature_type)
            if feature_type == 'source':
                strain = qualifiers.get('isolate')
                # print(strain)
                if strain:
                    strain_name_lst.append(strain[0])
                else:
                    strain_name_lst.append('unknown')

                collect_date = qualifiers.get('collection_date')
                if collect_date:
                    date_lst.append(collect_date[0])
                else:
                    date_lst.append('unknown')

                country = qualifiers.get('country')
                if country:
                    country_lst.append(country[0].split(':')[0])
                else:
                    country_lst.append('unknown')

                location = feature.location
                seq = location.extract(record.seq)
                seq_lst.append(seq)
                seq_len_lst.append(len(seq))

    df = pd.DataFrame()
    df['full_id'] = seq_id_lst
    df['strain_name'] = strain_name_lst
    df['seq'] = seq_lst
    df['seq_len'] = seq_len_lst
    df['date'] = date_lst
    df['country'] = country_lst
    df.to_csv(path + file_obj[:-3] + '.csv', index=False)

    WriteFastaAdd(path + file_obj[:-3] + '.fasta', seq_id_lst, seq_lst)
