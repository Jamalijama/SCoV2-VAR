import numpy as np

# import Levenshtein as ls
# coarse-grained read end mapping

'''
dnts: deca-nts, 10 nts
ents: ennea-nts, 9 nts 
onts: oct-nts, 8 nts
snts: sept-nts, 7 nts
hnts: hex-nts, 6 nts
pnts: penta-nts, 5 nts
'''
dict_xnt = {'dnt': 10, 'ent': 9, 'ont': 8, 'snt': 7, 'hnt': 6, 'pnt': 5}

lst_deg_base = ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
dict_deg_base = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'S': ['G', 'C'],
                 'W': ['A', 'T'], 'H': ['A', 'T', 'C'], 'B': ['G', 'T', 'C'], 'V': ['G', 'A', 'C'],
                 'D': ['G', 'A', 'T'], 'N': ['A', 'T', 'C', 'G']}


def mapping5_xnt(seq, xnt, dict_locs_xnts, skip0=0):
    # find the first start location in reference sequence
    flag, locsMapped5ter = False, []

    nts_num = dict_xnt[xnt]

    # Skip skip0 bases
    xnt_5ter0 = seq[skip0:skip0 + nts_num]
    xnt_5ter1 = seq[skip0 + nts_num:skip0 + nts_num * 2]

    if len(xnt_5ter0) != nts_num or len(xnt_5ter1) != nts_num:
        return flag, locsMapped5ter

    # # have degenerate base
    # for i in xnt_5ter0:
    #     if i in lst_deg_base:
    #         return flag, locsMapped5ter
    #
    # for i in xnt_5ter1:
    #     if i in lst_deg_base:
    #         return flag, locsMapped5ter

    # Get the positions of xnt_5ter0 and xnt_5ter1 in reference sequence
    locs_xnt_5ter0 = dict_locs_xnts[xnt_5ter0]
    locs_xnt_5ter1 = dict_locs_xnts[xnt_5ter1]

    # Get positional differences of xnt_5ter1 and xnt_5ter0 in reference sequence
    dlt_locs_lst = [loc_5ter1 - loc_5ter0
                    for loc_5ter0 in locs_xnt_5ter0
                    for loc_5ter1 in locs_xnt_5ter1]

    # Get the positions list of xnt_5ter0 in reference sequence
    loc_5ter0_lst = [loc_5ter0 for loc_5ter0 in locs_xnt_5ter0]

    for i in range(len(dlt_locs_lst)):

        # Get the i positional difference
        locdlt = dlt_locs_lst[i]

        # If the number of positional differences is equal to the number of bases in the small fragment
        if locdlt == nts_num:
            # matched
            flag = True
            start_loc = loc_5ter0_lst[int(i / len(locs_xnt_5ter1))]
            locsMapped5ter.append(start_loc)

            # Stop when find the first one.
            break

    return flag, locsMapped5ter


def mapping3_xnt(seq, xnt, loc_Mapped5ter, dict_locs_xnts, skip1=0):
    # find the first ending location in reference sequence
    # loc_Mapped5ter, an int, the location of the first mapped base in 5' terminal
    # locsMapped5ter, a list of ints, the locations of the mapped bases in 5' terminal

    nts_num = dict_xnt[xnt]

    readlen = len(seq)

    # the range of reference sequence where the aligning sequence may lie
    loc0_segref = loc_Mapped5ter + int(readlen * 0.8)
    loc1_segref = loc_Mapped5ter + int(readlen * 1.2)

    flagMapped3ter, locMapped3ter = False, 0

    # get the ending x bases segment and their start locations in reference sequence
    if skip1 == 0:
        xnt_3ter0 = seq[-(nts_num + skip1):]
    else:
        xnt_3ter0 = seq[-(nts_num + skip1):-skip1]
    xnt_3ter1 = seq[-(nts_num * 2 + skip1):-(nts_num + skip1)]

    if len(xnt_3ter0) != nts_num or len(xnt_3ter1) != nts_num:
        return flagMapped3ter, locMapped3ter

    locs_xnt_3ter0 = dict_locs_xnts[xnt_3ter0]
    locs_xnt_3ter1 = dict_locs_xnts[xnt_3ter1]

    # do not have such segment, not mapped
    if (len(locs_xnt_3ter0) != 0) and (len(locs_xnt_3ter1) != 0):
        # get the locations list in the range of reference sequence
        locs_xnt_3ter0 = [locc for locc in locs_xnt_3ter0
                          if (locc >= loc0_segref) & (locc <= loc1_segref)]
        locs_xnt_3ter1 = [locc for locc in locs_xnt_3ter1
                          if (locc >= loc0_segref) & (locc <= loc1_segref)]

        dlt_locs3ter_lst = [loc0 - loc1
                            for loc0 in locs_xnt_3ter0
                            for loc1 in locs_xnt_3ter1]
        loc1_xnt_3ter0_lst = [loc0 for loc0 in locs_xnt_3ter0]

        for j in range(len(dlt_locs3ter_lst)):
            dlt_locs3ter = dlt_locs3ter_lst[j]

            if dlt_locs3ter == nts_num:
                flagMapped3ter = True
                locMapped3ter = loc1_xnt_3ter0_lst[int(j / len(locs_xnt_3ter1))] + nts_num
                break

            else:
                flagMapped3ter = False
                locMapped3ter = 0

    return flagMapped3ter, locMapped3ter


# coarse-grained read end mapping
# fine-/coarse-grained read aligning


def ham_distance0(seq0, seq1):
    '''
    Parameters
    ----------
    seq0 : String
        DESCRIPTION:
            One of the two sequences to calculate Hamming distance.
    seq1 : String
        DESCRIPTION:
            The other of the two sequences to calculate Hamming distance.

    Returns
    -------
    TYPE: int
        DESCRIPTION:
            A Hamming editing distance of between two sequences.
    '''
    count = 0
    for xi, yi in zip(seq0, seq1):
        if xi != yi:
            if yi in lst_deg_base:
                lst_yi = dict_deg_base[yi]
                if xi not in lst_yi:
                    count += 1
            else:
                count += 1

    return count


def ham_distance1(seq0, seq1):
    # segmentation for comparison the hamming distance
    '''
        Parameters
        ----------
        seq0 : String
            DESCRIPTION:
                One of the two sequences to calculate Hamming distance.
        seq1 : String
            DESCRIPTION:
                The other of the two sequences to calculate Hamming distance.

        Returns
        -------
        TYPE: int
            DESCRIPTION:
                A Hamming editing distance of between two sequences.
    '''
    count = 0
    len_min_seq = min(len(seq0), len(seq1))
    for nt_i in range(0, len_min_seq - 20, 20):
        segseq0 = seq0[nt_i:nt_i + 20]
        segseq1 = seq1[nt_i:nt_i + 20]
        if segseq0 != segseq1:
            count0 = 0
            for xi, yi in zip(segseq0, segseq1):
                if xi != yi:
                    if yi in lst_deg_base:
                        lst_yi = dict_deg_base[yi]
                        if xi not in lst_yi:
                            count0 += 1
                    else:
                        count0 += 1
            count += count0
    return count


def lev_distance(seq0, seq1):
    '''
    Parameters
    ----------
    seq0 : String
        DESCRIPTION:
            One of the two sequences to calculate Levenshtein distance.
    seq1 : String
        DESCRIPTION:
            The other of the two sequences to calculate Levenshtein distance.

    Returns
    -------
    TYPE: int.
        DESCRIPTION:
            A Levenshtein distance of int type between two sequences.
    '''

    if len(seq0) >= len(seq1):
        seqL, seqS, seqlen = seq0, seq1, len(seq0)
    else:
        seqL, seqS, seqlen = seq1, seq0, len(seq1)

    count_arr = np.arange(seqlen + 1)

    for i in range(1, len(seqL) + 1):
        ref1 = count_arr[0]
        count_arr[0] += 1

        for j in range(1, len(seqS) + 1):
            ref2 = count_arr[j]

            if seqL[i - 1] == seqS[j - 1]:
                count_arr[j] = ref1
            else:
                count_arr[j] = min(ref1, min(count_arr[j - 1], count_arr[j])) + 1
            ref1 = ref2

    return count_arr[len(seqS)]


# nt_lst = ['-', 'T', 'C', 'A', 'G']
# dnt_lst = [i + j for i in nt_lst for j in nt_lst]
# dnt_score_lst = [
#     0, 0, 0, 0, 0,
#     0, 4, -1, -1, -1,
#     0, -1, 4, -1, -1,
#     0, -1, -1, 4, -1,
#     0, -1, -1, -1, 4
# ]
# dnt_score_dict = dict(zip(dnt_lst, dnt_score_lst))

nt_lst = ['-', 'T', 'C', 'A', 'G', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
dnt_lst = [i + j for i in nt_lst for j in nt_lst]
dnt_score_lst = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 4, -1, -1, -1, -1, 4, -1, 4, -1, 4, 4, 4, -1, 4, 4,
    0, -1, 4, -1, -1, -1, 4, 4, -1, 4, -1, 4, 4, 4, -1, 4,
    0, -1, -1, 4, -1, 4, -1, 4, -1, -1, 4, 4, -1, 4, 4, 4,
    0, -1, -1, -1, 4, 4, -1, -1, 4, 4, -1, -1, 4, 4, 4, 4,
    0, -1, -1, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 4, 4, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, -1, 4, 4, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 4, -1, -1, 4, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1,
    0, -1, 4, -1, 4, -1, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1,
    0, 4, -1, 4, -1, -1, -1, -1, -1, -1, 4, -1, -1, -1, -1, -1,
    0, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, 4, -1, -1, -1, -1,
    0, 4, 4, -1, 4, -1, -1, -1, -1, -1, -1, -1, 4, -1, -1, -1,
    0, -1, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, 4, -1, -1,
    0, 4, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4, -1,
    0, 4, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4,
]
dnt_score_dict = dict(zip(dnt_lst, dnt_score_lst))


def get_matrix_max(matrix, seq1, seq2):
    # get the index of max_score
    Max = matrix.max()

    for i in range(len(seq2), -1, -1):
        for j in range(len(seq1), -1, -1):
            if matrix[i][j] == Max:
                return i, j


def align_SW(seq1, seq2, dnt_dict=dnt_score_dict):
    # get aligned string of seq1 and seq2 based on local optimum solution
    '''
    Parameters
    ----------
    seq1 : String.
        Read/Ref sequence.
    seq2 : String.
        Ref/Read sequence.
    dnt_dict : Optional.
        DESCRIPTION. The default is dnt_score_dict.
        Dictionary of score of each 2-nt long "primer" for the ref seq,
        to align two seqs with the algorithm of Smith-Waterman.

    Returns
    -------
    aligned_seq1 / aligned_seq2: String.
        The read/ref seq post alignment.
    '''
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    s1, s2, gap = '', '', -2
    best_matrix = np.empty(shape=(len(seq2) + 1, len(seq1) + 1), dtype=int)

    # calculate the score matrix
    for i in range(len(seq2) + 1):
        for j in range(len(seq1) + 1):
            if (i == 0) or (j == 0):
                best_matrix[i][j] = 0
            else:
                match = dnt_dict[seq2[i - 1] + seq1[j - 1]]
                gap1_score = best_matrix[i - 1][j] + gap
                gap2_score = best_matrix[i][j - 1] + gap
                match_score = best_matrix[i - 1][j - 1] + match
                score = max(gap1_score, gap2_score, match_score)
                if score > 0:
                    best_matrix[i][j] = score
                else:
                    best_matrix[i][j] = 0

    # traceback to find the aligned sequence
    i, j = get_matrix_max(best_matrix, seq1, seq2)

    # Padding the end of the two sequences
    k = len(seq1)
    while k > j:
        s1 += seq1[k - 1]
        s2 += '-'
        k -= 1

    while best_matrix[i][j] != 0:
        match = dnt_dict[seq2[i - 1] + seq1[j - 1]]
        if i > 0 and j > 0 and best_matrix[i][j] == best_matrix[i - 1][j - 1] + match:
            s1 += seq1[j - 1]
            s2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif i > 0 and best_matrix[i, j] == best_matrix[i - 1, j] + gap:
            s1 += '-'
            s2 += seq2[i - 1]
            i -= 1
        else:
            s1 += seq1[j - 1]
            s2 += '-'
            j -= 1

    # Padding the start of the two sequences
    while j > 0:
        s1 += seq1[j - 1]
        s2 += '-'
        j -= 1

    aligned_seq1, aligned_seq2 = s1[::-1], s2[::-1]

    return aligned_seq1, aligned_seq2  # , best_matrix


seq1 = 'ACCACATTCGAAGTTTACCGTCGACCTCTCGCGATCTCCGAGCTACGAGACGTCGTACGTAGC'
seq2 = 'ACCACATTCACCGTCGACCTCTCGCGAACGTCGTACGTAGC'
s1, s2 = align_SW(seq1, seq2)


# print(seq1)
# print(s1)
# print(s2)


# fine-/coarse-grained read aligning
# multi-grained reference seq parse


def ref_dnt_parse(ref):
    # get a dictionary of locations for each of serially-cut dnt from a reference
    '''
    Parameters
    ----------
    ref : String.
        Reference sequence for alignment.

    Returns
    -------
    dict_locs_dnts (dnt: 10-nts): Dictionary
        A dictionary of locations for each of serially-cut dnt from a reference.
    '''

    seqlenTM = len(ref)

    nts = ['T', 'C', 'A', 'G']
    # generate the ten bases list
    dnt_lst = [nt0 + nt1 + nt2 + nt3 + nt4 + nt5 + nt6 + nt7 + nt8 + nt9
               for nt0 in nts for nt1 in nts for nt2 in nts
               for nt3 in nts for nt4 in nts for nt5 in nts
               for nt6 in nts for nt7 in nts for nt8 in nts
               for nt9 in nts]

    # a list contains list of 4**10 numbers
    locs_dnts_list = [[]] * (4 ** 10)

    dict_locs_dnts = dict(zip(dnt_lst, locs_dnts_list))

    for i in range(0, seqlenTM - 10 + 1, 1):
        dnt_i = ref[i: i + 10]
        lst_locs_dnt_i = dict_locs_dnts[dnt_i]
        lst_locs_dnt_i = lst_locs_dnt_i + [i]  # list append
        dict_locs_dnts[dnt_i] = lst_locs_dnt_i

    return dict_locs_dnts


def ref_pnt_parse(ref):
    '''
    Parameters
    ----------
    ref : String.
        Reference sequence for alignment.

    Returns
    -------
    dict_locs_pnts (pnt: 5-nts) : Dictionary
        A dictionary of locations for each of serially-cut pnt from a reference.
    '''

    seqlenTM = len(ref)
    nts = ['T', 'C', 'A', 'G']

    pnt_lst = [nt0 + nt1 + nt2 + nt3 + nt4
               for nt0 in nts for nt1 in nts for nt2 in nts
               for nt3 in nts for nt4 in nts]

    locs_pnts_list = [[]] * (4 ** 5)
    dict_locs_pnts = dict(zip(pnt_lst, locs_pnts_list))

    for i in range(0, seqlenTM - 5 + 1, 1):
        pnt_i = ref[i: i + 5]
        lst_locs_pnt_i = dict_locs_pnts[pnt_i]
        lst_locs_pnt_i = lst_locs_pnt_i + [i]
        dict_locs_pnts[pnt_i] = lst_locs_pnt_i
    return dict_locs_pnts


def ShortCutAlign(ref, seq, seqlen_5ter2map):
    snpflag = False
    inflag = False
    delflag = False

    read_mapaligned = ''
    ref_mapaligned = ''

    inLoc_lst = []
    inNts_lst = []
    cutSeq_ntslst = []

    seqlen = len(seq)
    reflen = len(ref)

    ham_score = ham_distance0(ref, seq)
    freq_ham_score = ham_score / seqlen

    # SNP
    if freq_ham_score <= 0.05:
        snpflag = True
        if reflen == seqlen:
            read_mapaligned = seq
            ref_mapaligned = ref
        elif reflen < seqlen:
            for i in range(reflen, seqlen):
                inLoc_lst.append(i)
                inNts_lst.append(seq[i])
            read_mapaligned = seq[:reflen]
            ref_mapaligned = ref
        else:
            for i in range(seqlen, reflen):
                seq += '-'
            read_mapaligned = seq
            ref_mapaligned = ref

    # another unspecified condition is that freq_HAM_score > 0.05
    elif (ham_distance0(ref[:seqlen_5ter2map],
                        seq[:seqlen_5ter2map]) < seqlen_5ter2map * 0.5):

        len_refseq = min(reflen, seqlen)

        # following scripts: to align read and ref_bimapped.

        interval = 9
        start_to_SW = -1
        for i in range(0, len_refseq - 6 + 1, 1):
            # the start num_xnt0 nts matched between ref&read
            segref0, segread00 = ref[i:i + 6], seq[i:i + 6]
            score_ham = ham_distance0(segref0, segread00)

            if score_ham >= 3:
                start_to_SW = i - interval if i - interval > 0 else 0
                break

        end_to_SW = 0
        for j in range(0, len_refseq - start_to_SW - 6, 1):
            # the last num_xnt1 nts matched between ref&read

            ref2_local3, read00_local3 = ref[-j - 6:-j], seq[-j - 6:-j]
            score_ham3 = ham_distance0(ref2_local3, read00_local3)

            if score_ham3 >= 3:
                end_to_SW = -(j - interval) if -(j - interval) < -1 else -1
                break

        if (start_to_SW != -1) & (end_to_SW != 0):
            start_to_SW -= int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
            end_to_SW += int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
            if start_to_SW < 0:
                start_to_SW = 0
            if end_to_SW > -1:
                end_to_SW = -1
            ref_for_SW = ref[start_to_SW: end_to_SW]
            seq_for_SW = seq[start_to_SW: end_to_SW]

            # Smith-Waterman for seq without 'N'
            # seq_for_SW = seq_for_SW.replace('N', '-')
            ref_aligned0, seq_aligned0 = align_SW(ref_for_SW, seq_for_SW)

            ref_mapped_aligned = ref[: start_to_SW] + ref_aligned0 + ref[end_to_SW:]
            seq_mapped_aligned0 = seq[: start_to_SW] + seq_aligned0 + seq[end_to_SW:]

            if seq_aligned0.count('-') > 0:
                delflag = True

            # only SNP
            if reflen == len(ref_mapped_aligned):
                read_mapaligned = seq_mapped_aligned0
                ref_mapaligned = ref_mapped_aligned
                inflag = False
            # insert
            elif reflen < len(ref_mapped_aligned):
                cutSeq_ref = []
                for nt_i in range(len(ref_mapped_aligned)):
                    nt_ref = ref_mapped_aligned[nt_i]
                    nt_seq = seq_mapped_aligned0[nt_i]

                    if nt_ref == '-' and nt_seq == '-':
                        continue
                    elif nt_ref == '-':
                        inLoc_lst.append(nt_i)
                        inNts_lst.append(nt_seq)
                    else:
                        cutSeq_ref.append(nt_ref)
                        cutSeq_ntslst.append(nt_seq)

                read_mapaligned = ''.join(cutSeq_ntslst)
                ref_mapaligned = ''.join(cutSeq_ref)
                inflag = True
            # delete
            else:
                ref_mapaligned = ref_mapped_aligned
                read_mapaligned = seq_mapped_aligned0
                for _ in range(reflen - len(ref_mapped_aligned)):
                    read_mapaligned += '-'

    return snpflag, inflag, delflag, read_mapaligned, ref_mapaligned, inLoc_lst, inNts_lst


# multi-grained reference seq parse
# main FastAlign


def FastAlign(seq, ref, xnt0, xnt1, dict_locs_xnts0, dict_locs_xnts1, seqlen_5ter2map, segLen):
    read_mapaligned = ''
    ref_mapaligned = ''

    seqlen = len(seq)
    # print(seqlen)
    inLoc_lst, inNts_lst = [], []
    start_loc = -1
    end_loc = 0

    num_xnt0 = dict_xnt[xnt0]
    num_xnt1 = dict_xnt[xnt1]

    ham_score = seqlen
    # default is the maximum unmatched value.

    mloc = 0
    # segLen = 500
    mapped_start_loc_lst = []
    skip_lst = []

    for mloc in range(0, seqlen, segLen):
        for skip0 in range(0, 91, 5):
            ss = seq[mloc + skip0:mloc + skip0 + num_xnt0 * 2]
            flag = 0
            for i in range(len(ss)):
                if ss[i] in lst_deg_base:
                    flag = 1
            if flag == 1:
                continue
            else:
                mapping50 = mapping5_xnt(seq[mloc:], xnt0, dict_locs_xnts0, skip0=skip0)

                if mapping50[0]:  # mapping50 return flagMapped5ter, locsMapped5ter
                    skip_lst.append(skip0)
                    start_loc_lst0 = mapping50[1]

                    # following to filter start_loc
                    if len(start_loc_lst0) == 1:
                        mapped_start_loc_lst.append(start_loc_lst0[0])

                    elif len(start_loc_lst0) > 1:
                        ham_score_lst0 = []
                        for start_loc0 in start_loc_lst0:
                            ref1_local = ref[start_loc0: start_loc0 + seqlen - mloc]
                            ham_score0 = ham_distance1(ref1_local, seq[mloc + skip0:])
                            ham_score_lst0.append(ham_score0)

                        ham_score = min(ham_score_lst0)
                        index_ham_score = ham_score_lst0.index(ham_score)
                        mapped_start_loc_lst.append(start_loc_lst0[index_ham_score])

                    if len(mapped_start_loc_lst) >= 2:
                        ref1 = ref[mapped_start_loc_lst[-2]:mapped_start_loc_lst[-1] + num_xnt0]
                        seq1 = seq[mloc - segLen + skip_lst[-2]:mloc + skip_lst[-1] + num_xnt0]

                        ham_score = ham_distance0(ref1, seq1)
                        freq_ham_score = ham_score / len(seq1)

                        # only have SNP
                        if freq_ham_score <= 0.05 and len(ref1) == len(seq1):  # an empiric threshold, parameterizable
                            read_mapaligned += seq1[:-num_xnt0]
                            ref_mapaligned += ref1[:-num_xnt0]

                        # not only have snp
                        elif (ham_distance0(ref1[:seqlen_5ter2map],
                                            seq1[:seqlen_5ter2map]) < seqlen_5ter2map):
                            ref1len = len(ref1)
                            seq1len = len(seq1)
                            len_refseq = min(ref1len, seq1len)

                            interval = 9
                            start_to_SW = -1
                            for i in range(num_xnt0, len_refseq - 6 + 1, 1):
                                # the first num_xnt0 nts matched between ref and read

                                ref2_local5, read00_local5 = ref1[i:i + 6], seq1[i:i + 6]
                                score_ham5 = ham_distance0(ref2_local5, read00_local5)

                                # location i has insert or delete
                                if score_ham5 >= 3:  # an empirac threshold.
                                    start_to_SW = i - interval if i - interval > 0 else 0
                                    # start location for SW move forward 9 bases
                                    break

                            end_to_SW = 0
                            for j in range(num_xnt0, len_refseq - start_to_SW - 6, 1):
                                # the last num_xnt1 nts matched between ref and read
                                ref2_local3, read00_local3 = ref1[-j - 6:-j], seq1[-j - 6:-j]
                                score_ham3 = ham_distance0(ref2_local3, read00_local3)

                                if score_ham3 >= 3:
                                    end_to_SW = -(j - interval) if -(j - interval) < -1 else -1
                                    # ending location for SW move back 9 bases
                                    break

                            if (start_to_SW != -1) & (end_to_SW != 0):
                                start_to_SW -= int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
                                end_to_SW += int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
                                if start_to_SW < 0:
                                    start_to_SW = 0
                                if end_to_SW > -1:
                                    end_to_SW = -1
                                ref_for_SW = ref1[start_to_SW: end_to_SW]
                                seq_for_SW = seq1[start_to_SW: end_to_SW]
                                # length of ref_for_SW maybe not same as seq_for_SW
                                # seq_for_SW = seq_for_SW.replace('N', '-')  # Smith-Waterman for seq without 'N'
                                ref_aligned0, seq_aligned0 = align_SW(ref_for_SW, seq_for_SW)

                                ref_mapped_aligned0 = ref1[: start_to_SW] + ref_aligned0 + ref1[end_to_SW:]
                                seq_mapped_aligned0 = seq1[: start_to_SW] + seq_aligned0 + seq1[end_to_SW:]
                                len_ref_mapaligned = len(ref_mapped_aligned0)
                                len_ref1 = len(ref1)

                                # do not have insert
                                if len_ref1 == len_ref_mapaligned:
                                    read_mapaligned += seq_mapped_aligned0[:-num_xnt0]
                                    ref_mapaligned += ref_mapped_aligned0[:-num_xnt0]

                                # have insert
                                elif len_ref1 < len_ref_mapaligned:
                                    cutSeq_ntslst = []
                                    cutSeq_ref = []
                                    for nt_i in range(len_ref_mapaligned):
                                        nt_ref = ref_mapped_aligned0[nt_i]
                                        nt_seq = seq_mapped_aligned0[nt_i]

                                        if nt_ref == '-' and nt_seq == '-':
                                            continue
                                        elif nt_ref == '-':
                                            inLoc_lst.append(nt_i + mapped_start_loc_lst[-2])
                                            inNts_lst.append(nt_seq)
                                        else:
                                            cutSeq_ref.append(nt_ref)
                                            cutSeq_ntslst.append(nt_seq)

                                    read_mapaligned += ''.join(cutSeq_ntslst)[:-num_xnt0]
                                    ref_mapaligned += ''.join(cutSeq_ref)[:-num_xnt0]

                    # mapped, turn to next cut
                    break

    start_loc = mapped_start_loc_lst[0]
    skip0 = skip_lst[0]
    for skip1 in range(0, 61, 5):
        # following to filter end_loc
        if skip1 == 0:
            ss = seq[-(num_xnt1 * 2 + skip1):]
        else:
            ss = seq[-(num_xnt1 * 2 + skip1):-skip1]
        flag = 0
        for i in range(num_xnt1 * 2):
            if ss[i] in lst_deg_base:
                flag = 1
        if flag == 1:
            continue
        else:
            mapping30 = mapping3_xnt(seq[mloc:], xnt1, mapped_start_loc_lst[-1], dict_locs_xnts1, skip1=skip1)

            # following to map the read 3' terminal
            if mapping30[0]:
                end_loc = mapping30[1]

                if end_loc != 0:
                    ref1 = ref[mapped_start_loc_lst[-1]:end_loc]
                    if skip1 == 0:
                        seq1 = seq[mloc + skip_lst[-1]:]
                    else:
                        seq1 = seq[mloc + skip_lst[-1]:-skip1]
                    # above to filter end_loc
                    ham_score = ham_distance0(ref1, seq1)
                    freq_ham_score = ham_score / len(seq1)

                    # only have SNP
                    if freq_ham_score <= 0.05 and len(ref1) == len(seq1):  # an empiric threshold, parameterizable
                        read_mapaligned += seq1
                        ref_mapaligned += ref1

                        read_mapaligned1 = ''
                        ref_mapaligned1 = ''
                        if start_loc > 0 and skip0 != 0:
                            mappingbegin = ShortCutAlign(ref[:start_loc][::-1], seq[:skip0][::-1],
                                                         seqlen_5ter2map)
                            if len(mappingbegin[5]) != 0:
                                for i in range(len(mappingbegin[5]) - 1, -1, -1):
                                    inLoc_lst.append(int(start_loc - mappingbegin[5][i] - 1))
                                    inNts_lst.append(mappingbegin[6][i])

                            # move forward the start location
                            start_loc -= len(mappingbegin[4])

                            read_mapaligned1 = mappingbegin[3][:len(mappingbegin[4])][
                                               ::-1] + read_mapaligned
                            ref_mapaligned1 = mappingbegin[4][::-1] + ref_mapaligned

                        if read_mapaligned1 != '':
                            read_mapaligned = read_mapaligned1
                            ref_mapaligned = ref_mapaligned1

                        read_mapaligned1 = ''
                        ref_mapaligned1 = ''
                        if end_loc < len(ref) and skip1 != 0:
                            mappingend = ShortCutAlign(ref[end_loc:], seq[-skip1:], seqlen_5ter2map)

                            if len(mappingend[5]) != 0:
                                for i in range(len(mappingend[5])):
                                    inLoc_lst.append(int(end_loc + mappingend[5][i] + 1))
                                    inNts_lst.append(mappingend[6][i])

                            # move back the ending location
                            end_loc += len(mappingend[4])

                            read_mapaligned1 = read_mapaligned + mappingend[3][:len(mappingend[4])]
                            ref_mapaligned1 = ref_mapaligned + mappingend[4]

                        if read_mapaligned1 != '':
                            read_mapaligned = read_mapaligned1
                            ref_mapaligned = ref_mapaligned1

                        if start_loc > 0:
                            for i in range(start_loc):
                                read_mapaligned = '-' + read_mapaligned
                                ref_mapaligned = ref[start_loc - i - 1] + ref_mapaligned

                        if end_loc < len(ref):
                            for i in range(end_loc, len(ref)):
                                read_mapaligned += '-'
                                ref_mapaligned += ref[i]

                        return read_mapaligned, start_loc, end_loc, inLoc_lst, inNts_lst, ref_mapaligned

                    # not only have snp
                    elif (ham_distance0(ref1[:seqlen_5ter2map],
                                        seq1[:seqlen_5ter2map]) < seqlen_5ter2map * 0.5):
                        ref1len = len(ref1)
                        seq1len = len(seq1)
                        len_refseq = min(ref1len, seq1len)

                        interval = 9
                        start_to_SW = -1
                        for i in range(num_xnt0, len_refseq - 6 + 1, 1):
                            # the first num_xnt0 nts matched between ref and read

                            ref2_local5, read00_local5 = ref1[i:i + 6], seq1[i:i + 6]
                            score_ham5 = ham_distance0(ref2_local5, read00_local5)

                            # location i has insert or delete
                            if score_ham5 >= 3:  # an empirac threshold.
                                # start_t = i + 6
                                start_to_SW = i - interval if i - interval > 0 else 0
                                # start location for SW move forward 9 bases
                                break

                        end_to_SW = 0
                        # end_t = 0
                        for j in range(num_xnt1, len_refseq - start_to_SW - 6, 1):
                            # the last num_xnt1 nts matched between ref and read

                            ref2_local3, read00_local3 = ref1[-j - 6:-j], seq1[-j - 6:-j]
                            score_ham3 = ham_distance0(ref2_local3, read00_local3)

                            if score_ham3 >= 3:
                                end_to_SW = -(j - interval) if -(j - interval) < -1 else -1
                                # ending location for SW move back 9 bases
                                break

                        if (start_to_SW != -1) & (end_to_SW != 0):
                            start_to_SW -= int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
                            end_to_SW += int((len(seq) - start_to_SW + end_to_SW - 2 * interval) / interval)
                            if start_to_SW < 0:
                                start_to_SW = 0
                            if end_to_SW > -1:
                                end_to_SW = -1
                            ref_for_SW = ref1[start_to_SW: end_to_SW]
                            seq_for_SW = seq1[start_to_SW: end_to_SW]
                            # length of ref_for_SW maybe not same as seq_for_SW
                            # seq_for_SW = seq_for_SW.replace('N', '-')  # Smith-Waterman for seq without 'N'
                            ref_aligned0, seq_aligned0 = align_SW(ref_for_SW, seq_for_SW)

                            ref_mapped_aligned0 = ref1[: start_to_SW] + ref_aligned0 + ref1[end_to_SW:]
                            seq_mapped_aligned0 = seq1[: start_to_SW] + seq_aligned0 + seq1[end_to_SW:]
                            len_ref_mapaligned = len(ref_mapped_aligned0)
                            len_ref1 = len(ref1)

                            # do not have insert
                            if len_ref1 == len_ref_mapaligned:
                                read_mapaligned += seq_mapped_aligned0
                                ref_mapaligned += ref_mapped_aligned0

                            # have insert
                            elif len_ref1 < len_ref_mapaligned:
                                cutSeq_ntslst = []
                                cutSeq_ref = []
                                for nt_i in range(len_ref_mapaligned):
                                    nt_ref = ref_mapped_aligned0[nt_i]
                                    nt_seq = seq_mapped_aligned0[nt_i]

                                    if nt_ref == '-' and nt_seq == '-':
                                        continue
                                    elif nt_ref == '-':
                                        inLoc_lst.append(nt_i + mapped_start_loc_lst[-1] + len(inNts_lst))
                                        inNts_lst.append(nt_seq)
                                    else:
                                        cutSeq_ref.append(nt_ref)
                                        cutSeq_ntslst.append(nt_seq)

                                read_mapaligned += ''.join(cutSeq_ntslst)
                                ref_mapaligned += ''.join(cutSeq_ref)

                            read_mapaligned1 = ''
                            ref_mapaligned1 = ''
                            if start_loc > 0 and skip0 != 0:
                                mappingbegin = ShortCutAlign(ref[:start_loc][::-1], seq[:skip0][::-1],
                                                             seqlen_5ter2map)
                                if len(mappingbegin[5]) != 0:
                                    for i in range(len(mappingbegin[5]) - 1, -1, -1):
                                        inLoc_lst.append(int(start_loc - mappingbegin[5][i] - 1))
                                        inNts_lst.append(mappingbegin[6][i])

                                # move forward the start location
                                start_loc -= len(mappingbegin[4])

                                read_mapaligned1 = mappingbegin[3][:len(mappingbegin[4])][
                                                   ::-1] + read_mapaligned
                                ref_mapaligned1 = mappingbegin[4][::-1] + ref_mapaligned

                            if read_mapaligned1 != '':
                                read_mapaligned = read_mapaligned1
                                ref_mapaligned = ref_mapaligned1

                            read_mapaligned1 = ''
                            ref_mapaligned1 = ''
                            if end_loc < len(ref) and skip1 != 0:
                                mappingend = ShortCutAlign(ref[end_loc:], seq[-skip1:], seqlen_5ter2map)

                                if len(mappingend[5]) != 0:
                                    for i in range(len(mappingend[5])):
                                        inLoc_lst.append(int(end_loc + mappingend[5][i] + 1))
                                        inNts_lst.append(mappingend[6][i])

                                # move back the ending location
                                end_loc += len(mappingend[4])

                                read_mapaligned1 = read_mapaligned + mappingend[3][:len(mappingend[4])]
                                ref_mapaligned1 = ref_mapaligned + mappingend[4]

                            if read_mapaligned1 != '':
                                read_mapaligned = read_mapaligned1
                                ref_mapaligned = ref_mapaligned1

                            if start_loc > 0:
                                for i in range(start_loc):
                                    read_mapaligned = '-' + read_mapaligned
                                    ref_mapaligned = ref[start_loc - i - 1] + ref_mapaligned

                            if end_loc < len(ref):
                                for i in range(end_loc, len(ref)):
                                    read_mapaligned += '-'
                                    ref_mapaligned += ref[i]

                            return read_mapaligned, start_loc, end_loc, inLoc_lst, inNts_lst, ref_mapaligned

    return read_mapaligned, start_loc, end_loc, inLoc_lst, inNts_lst, ref_mapaligned
