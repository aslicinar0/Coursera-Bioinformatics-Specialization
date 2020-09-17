import numpy as np

with open("BLOSUM62.txt","r") as blosum_file:
    x = blosum_file.readline().rstrip()
    aas = list(x[3:].split("  "))
    blosum_62 = np.matrix(np.loadtxt("BLOSUM62.txt", skiprows=1, usecols=range(1, 21)))


with open("week3-1-Alignment_with_Affine_Gap_Penalties.txt","r") as x:
    sequences= x.readlines()
    sequence_1 = sequences[0].strip()
    sequence_2 = sequences[1].strip()


def lcs_backtrack(v:str,w:str,blosum62: np.ndarray,aas:list):
    gap_op = 11 #gap opening penalty
    gap_ex = 1  #gap extention penalty
    len_v = len(v)
    len_w = len(w)
    s_lower = np.zeros((len_v + 1,len_w + 1))
    s_middle = np.zeros((len_v + 1, len_w + 1))
    s_upper = np.zeros((len_v + 1, len_w + 1))
    backtrack_middle = np.zeros((len_v + 1, len_w + 1))
    backtrack_lower = np.zeros((len_v + 1, len_w + 1))
    backtrack_upper = np.zeros((len_v + 1, len_w + 1))

    for i in range(1, len_v + 1):
        s_lower[i][0] = s_lower[0][0] - (gap_op + (i-1)*gap_ex)
        s_middle[i][0] = s_middle[0][0] - (gap_op + (i-1)*gap_ex)
        s_upper[i][0] = np.NINF

    for i in range(1, len_w + 1):
        s_upper[0][i] = s_upper[0][0] - (gap_op + (i - 1) * gap_ex)
        s_middle[0][i] = s_middle[0][0] - (gap_op + (i - 1) * gap_ex)
        s_lower[0][i] = np.NINF

    for i in range(1,len_v+1):
        for j in range(1, len_w+1):
            v_index = aas.index(v[i - 1])
            w_index = aas.index(w[j - 1])
            """
            if v[i - 1] == w[j - 1]:
                match_mismatch_point = 2
            else:
                match_mismatch_point = -1
            """
            match_mismatch_point = blosum62[v_index, w_index]

            s_lower[i][j] = max(s_lower[i-1][j]-gap_ex, s_middle[i-1][j]-gap_op)
            if s_lower[i][j] == s_lower[i-1][j]-gap_ex:
                backtrack_lower[i][j] = "2"
            else:
                backtrack_lower[i][j] = "5"

            s_upper[i][j] = max(s_upper[i][j-1]-gap_ex, s_middle[i][j-1]-gap_op)
            if s_upper[i][j] == s_upper[i][j-1]-gap_ex:
                backtrack_upper[i][j] = "6"
            else:
                backtrack_upper[i][j] = "5"
            s_middle[i][j] = max(s_lower[i][j], s_upper[i][j], s_middle[i-1][j-1] + match_mismatch_point)
            if s_middle[i][j] == s_lower[i][j]:
                backtrack_middle[i][j] = "2"
            elif s_middle[i][j] == s_upper[i][j]:
                backtrack_middle[i][j] = "6"
            else:
                backtrack_middle[i][j] = "5"
    backtrack_middle = np.delete(backtrack_middle, 0, 0)
    backtrack_middle = np.delete(backtrack_middle, 0, 1)
    backtrack_lower = np.delete(backtrack_lower, 0, 0)
    backtrack_lower = np.delete(backtrack_lower, 0, 1)
    backtrack_upper = np.delete(backtrack_upper, 0, 0)
    backtrack_upper = np.delete(backtrack_upper, 0, 1)


    print(s_middle[len_v, len_w])

    return backtrack_middle, backtrack_lower, backtrack_upper


backtracks = lcs_backtrack(sequence_2,sequence_1,blosum_62,aas)
backtrack_middle_2 = backtracks[0]
backtrack_lower_2 = backtracks[1]
backtrack_upper_2 = backtracks[2]


def b_middle(i:int, j: int, seq_1:str, seq_2:str):

    if i >= 0 and j >= 0:
        seq_1 += sequence_1[j]
        seq_2 += sequence_2[i]
        i -= 1
        j -= 1
        if backtrack_middle_2[i][j] == 5.0:
            b_middle(i,j,seq_1,seq_2)
        elif backtrack_middle_2[i][j] == 2.0:
            b_lower(i,j,seq_1,seq_2)
        elif backtrack_middle_2[i][j] == 6.0:
            b_upper(i,j,seq_1,seq_2)
    else:
        seq_1_res = ''.join(reversed(seq_1))
        seq_2_res = ''.join(reversed(seq_2))
        print(seq_1_res, seq_2_res)
        return


def b_lower(i :int, j: int, seq_1:str, seq_2:str):
    if i >= 0:
        if backtrack_lower_2[i][j] == 2.0:
            seq_1 += "-"
            seq_2 += sequence_2[i]
            i -= 1
            b_lower(i,j,seq_1,seq_2)
        else:
            seq_1 += "-"
            seq_2 += sequence_2[i]
            i -= 1
            b_middle(i,j,seq_1,seq_2)

    else:
        seq_1 += ''.join(reversed(sequence_1[:j+1]))
        seq_2 += (j+1) * "-"
        seq_1_res = ''.join(reversed(seq_1))
        seq_2_res = ''.join(reversed(seq_2))
        print(seq_1_res, seq_2_res)
        return

def b_upper(i:int, j: int, seq_1:str, seq_2:str):
    if j >= 0:
        if backtrack_upper_2[i][j] == 6.0:
            seq_1 += sequence_1[j]
            seq_2 += "-"
            j -= 1
            b_upper(i,j,seq_1,seq_2)
        else:
            seq_1 += sequence_1[j]
            seq_2 += "-"
            j -= 1
            b_middle(i,j,seq_1,seq_2)
    else:
        seq_1 += (i + 1) * "-"
        seq_2 += ''.join(reversed(sequence_2[:i + 1]))
        seq_1_res = ''.join(reversed(seq_1))
        seq_2_res = ''.join(reversed(seq_2))
        print(seq_1_res,seq_2_res)
        return


def iterative_output_lcs(v:str,w:str):

    i = len(v) - 1
    j = len(w) - 1
    seq_1 = ""
    seq_2 = ""

    if backtrack_middle_2[i][j] == 5:
        b_middle(i,j,seq_1,seq_2)
    elif backtrack_middle_2[i][j] == 2:
        b_lower(i,j,seq_1,seq_2)
    else:
        b_upper(i,j,seq_1,seq_2)

    return



iterative_output_lcs(sequence_2,sequence_1)

