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
    gap_op = 15 #gap opening penalty
    gap_ex = 5  #gap extention penalty
    len_v = len(v)
    len_w = len(w)
    s_lower = np.zeros((len_v + 1,len_w + 1))
    s_middle = np.zeros((len_v + 1, len_w + 1))
    s_upper = np.zeros((len_v + 1, len_w + 1))
    backtrack = np.zeros((len_v + 1, len_w + 1))

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
            if v[i - 1] == w[j - 1]:
                match_mismatch_point = 5
            else:
                match_mismatch_point = -2
            #match_mismatch_point = blosum62[v_index, w_index]
            s_lower[i][j] = max(s_lower[i-1][j]-gap_ex, s_middle[i-1][j]-gap_op)
            s_upper[i][j] = max(s_upper[i][j-1]-gap_ex, s_middle[i][j-1]-gap_op)
            s_middle[i][j] = max(s_lower[i][j], s_upper[i][j], s_middle[i-1][j-1] + match_mismatch_point)
            if s_middle[i][j] == s_lower[i][j]:
                backtrack[i][j] = "2"
            elif s_middle[i][j] == s_upper[i][j]:
                backtrack[i][j] = "6"
            else:
                backtrack[i][j] = "5"
    backtrack = np.delete(backtrack, 0, 0)
    backtrack = np.delete(backtrack, 0, 1)
    print(s_lower)
    print(s_middle)
    print(s_upper)
    print(s_middle[len_v, len_w])
    print(backtrack)
    return backtrack


def iterative_output_lcs(backtrack:np.ndarray, v:str,w:str):
    """
    a = 6*np.ones((1,len(w)))
    bk = np.vstack((a,backtrack))
    print(bk)
    b = 2*np.ones((len(v)+1,1))
    print(b)
    bk2=np.column_stack((b,bk))
    print(bk2)
    backtrack = bk2
    """
    lcs = ""
    i = len(v) - 1
    j = len(w) - 1
    seq_1 = ""
    seq_2 = ""
    while i >= 0 and j >= 0:
        if backtrack[i][j] == 5:
            i -= 1
            j -= 1
            seq_1 += w[j + 1]
            seq_2 += v[i + 1]
        elif backtrack[i][j] == 2:
            i -= 1
            seq_1 += "-"
            seq_2 += v[i + 1]
        else:
            j -= 1
            seq_1 += w[j + 1]
            seq_2 += "-"
    print(i,j)
    if i < 0 :
        seq_1 += ''.join(reversed(w[:j+1]))
        seq_2 += (j+1) * "-"
    if j < 0 :
        seq_1 += (i+1) * "-"
        seq_2 += ''.join(reversed(v[:i+1]))
    seq_1_res = ''.join(reversed(seq_1))
    seq_2_res = ''.join(reversed(seq_2))
    return seq_1_res, seq_2_res

print(iterative_output_lcs(lcs_backtrack(sequence_1,sequence_2,blosum_62,aas),sequence_1,sequence_2))

