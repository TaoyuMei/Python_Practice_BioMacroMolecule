from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

##  Part A: Define All functions needed
#   1. extract the base pairs for each sequence (2 functions)

def colle_mask_pair(seq_dict):
    """collect the index of every base pair from the pairing mask of a fasta file
    the input is a dictionary created by seq_dict = SeqIO.to_dict(seq_iterator)
    in the returned dictionary, each base pair appears twice as a: b and b: a respectively"""
    index_revindex = {}
    # key is the index of a base in the Pairingmask, value is the index of its pairing base in the reverse mask
    mask = seq_dict['Pairingmask'].seq
    rev_mask = seq_dict['Pairingmask'].seq[-1::-1]
    mask_len = len(mask)
    col_dif_group = {}  # a different number or character represents a different group of pairs (helix) in the mask
    for c in mask:  # collect all unique helix symbols
        if(c != '-'):  # dash means gaps, not matches
            col_dif_group[c] = c
    pair_num_list = list(col_dif_group.keys())
    print("all unique helix symbols", pair_num_list)

    for character in pair_num_list:
        begin_index = 0
        next_begin = 0
        for i in range(mask_len):
            if (mask[i] == character) and (begin_index < mask_len):
                for j in range(begin_index, mask_len):
                    if(mask[i] == rev_mask[j]):
                        index_revindex[i] = j
                        next_begin = j + 1
                        break
            begin_index = next_begin
#    print(mask[142], "match", rev_mask[index_revindex[142]])
    for key in index_revindex.keys():  # convert the index of the reverse mask to the index of the mask
        index_revindex[key] = mask_len - index_revindex[key] - 1
    print("base pair positions in the pairing mask: ", index_revindex)
    return index_revindex


def extract_pair(record, index_revindex):
    """extract all base pairs from a single RNA sequence
    the input is a value of a dictionary created by SeqIO.to_dict()
    and a dictionary containing index of matches in the Pairingmask;
    in the returned tuple list, each base pair appear only once"""
    pairs = {}  # positions and corresponding matching positions
    sequ = record.seq  # sequence of the record
    for index in index_revindex.keys():
        if(sequ[index].isupper() and sequ[index_revindex[index]].isupper()):
            pairs[index] = index_revindex[index]
    pairs_tuples_list = list(pairs.items())
    tmp_list = pairs_tuples_list
    for tupl in pairs_tuples_list:  # remove redundancy
        for tu in tmp_list:
            if(tupl[0] == tu[1] and tupl[1] == tu[0]):
                tmp_list.remove(tu)
                break
    pairs_tuples_list = tmp_list
    return pairs_tuples_list


#   2. calculate the base pair distance between every pairs of sequences and draw a histogram (1 function)

def cal_bp_distance(seq1: list, seq2: list):
    """calculate the base pair distance between the 2 given sequences
    the input is 2 lists of base pair index tuples"""
    bp_dist = 0
    # number of base pairs that are present in the sequence 1 but not in the sequence 2
    for bp in seq1:
        if bp not in seq2:
            bp_dist += 1
    # number of base pairs that are present in the sequence 2 but not in the sequence 1
    for bbp in seq2:
        if bbp not in seq1:
            bp_dist += 1
    return bp_dist


#   3. calculate the Hamming distance between every pairs of sequences and draw a histogram (1 function)

def cal_ham_dist(seq1: str, seq2: str):
    """calculate the Hamming distance between the 2 given sequences.
    the 2 sequences should have been aligned"""
    hmd = 0
    for m in range(len(seq1)):
        if (seq1[m] in ('a', 'u', 'g', 'c', 'A', 'U', 'G', 'C')):
            if(seq1[m].upper() != seq2[m].upper() and seq2[m] != 'n'):
                hmd += 1
        elif (seq1[m] == '-'):
            if(seq1[m] != seq2[m]):
                hmd += 1
        elif(seq2[m] == '-'):  # seq1[m] must be 'n'
            hmd += 1
    return hmd


## Part B: Main Part

if __name__ == "__main__":
    #   1. extract the base pairs for each sequence

    seq_iterator = SeqIO.parse("baci_tmrna.fasta", "fasta")  # parse the fasta file
    seq_dict = SeqIO.to_dict(seq_iterator)  # store every record in a dictionary, dict_name[seq_id].seq is the sequence
    print("the pairing mask: ", seq_dict['Pairingmask'].seq)
    print("the reserve pairing mask", seq_dict['Pairingmask'].seq[-1::-1])

    index_revindex = colle_mask_pair(seq_dict)

    record_pair_tuples = {}
    for record in seq_dict.values():  # get the base pairs in each sequences and store them in a dictionary
        if record.id != 'Pairingmask':
            record_pair_tuples[record.id] = extract_pair(record, index_revindex)
            print("sequence ID: ", record.id)
            print("sequence: ", record.seq)
            print("base pair positions: ", record_pair_tuples[record.id])

    #   2. calculate the base pair distance between every pairs of sequences and draw a histogram

    #  compute the base pair distance matrix
    print("The base pair distance matrix of the 12 tmRNA sequences: ")
    all_id = list(record_pair_tuples.keys())
    for id1 in all_id:
        seq_bpd = []  # collect the base pair distances between this sequence and any other sequences
        for id2 in all_id:
            bpd = cal_bp_distance(record_pair_tuples[id1], record_pair_tuples[id2])
            seq_bpd.append(bpd)
        bpd_array = np.array(seq_bpd)
        print(bpd_array)

    #  compute the pairwise base pair distances of all sequences and draw a histogram
    bpd_list = []
    all_id = list(record_pair_tuples.keys())
    for n in range(len(all_id) - 1):
        id1 = all_id.pop()
        for id2 in all_id:
            bpd = cal_bp_distance(record_pair_tuples[id1], record_pair_tuples[id2])
            bpd_list.append(bpd)
    plt.cla()
    plt.hist(np.array(bpd_list), facecolor="blue", edgecolor="black", density=1, bins=35)
    plt.title("the Distribution of Pairwise Base Pair Distance")
    plt.savefig('./the Distribution of Base Pair Distance.png', format='png')
    print("Base Pair Distances :", bpd_list)

    #   3. calculate the Hamming distance between every pairs of sequences and draw a histogram

    #  compute the Hamming distance matrix
    print("The Hamming distance matrix of the 12 tmRNA sequences: ")
    all_id = list(seq_dict.keys())
    all_id.pop(all_id.index("Pairingmask"))
    for id1 in all_id:
        seq_hmd = []  # collect the Hamming distances between this sequence and any other sequences
        for id2 in all_id:
            hmd = cal_ham_dist(seq_dict[id1].seq, seq_dict[id2].seq)
            seq_hmd.append(hmd)
        hmd_array = np.array(seq_hmd)
        print(hmd_array)

    #  compute the pairwise Hamming distances of all sequences and draw a histogram
    hmd_list = []
    all_id = list(seq_dict.keys())
    all_id.pop(all_id.index("Pairingmask"))
    for n in range(len(all_id) - 1):
        id1 = all_id.pop()
        for id2 in all_id:
            hmd = cal_ham_dist(seq_dict[id1].seq, seq_dict[id2].seq)
            hmd_list.append(hmd)
    plt.cla()
    plt.hist(np.array(hmd_list), facecolor="blue", density=1, bins=20, edgecolor="black")
    plt.title("the Distribution of Hamming Distance")
    plt.savefig('./the Distribution of Hamming Distance.png', format='png')
    print("Hamming Distances: ", hmd_list)


