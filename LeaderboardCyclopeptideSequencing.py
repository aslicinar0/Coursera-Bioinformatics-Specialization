"""
04.04.2020
Written by Aslı Gizem Çınar

This code finds the best sequence that can most closely create the given experimental spectrum(experimental_spec)
(which is a mass spectrometry outcome of an unknown cyclopeptide.) by using branch and bound algorithm.
(Should do this less than 3 minutes for a 22 aa long peptide.)
The branching step adds every amino acid to a peptide.(function starts with an empty peptide.)(accepts the aas with same mass as one. I and L, K and Q)
Bounding step trims the outcome of the branching step with respect to N.
"""

from operator import itemgetter
from itertools import product
import timeit

amino_acids =[]
aas_and_masses ={}
#aa_mass_for_cyclopeptide accepts the aas with same mass as one. I and L(showned as I), K and Q(showned as K)
#reads the txt file and creates a dict.
with open("aa_mass_for_cyclopeptide.txt") as fd:
    for line in fd.readlines():
        aa, sep, mass = line.strip().partition(" ")
        aas_and_masses[aa] = int(mass)
        amino_acids.append(aa)

# converts the given peptide sequence to its masses and returns list of them.
def aas_to_mass(peptide:str):
    aa_mass_list = []
    for aa in peptide:
        aa_mass_list.append(aas_and_masses[aa])
    return aa_mass_list

# creates the subspectrums of a given peptide in cyclic format.
def subspectrum_generator(peptide_str:str):
    peptide = aas_to_mass(peptide_str)
    L = len(peptide)
    looped = peptide + peptide[:L-2]
    return [0, sum(peptide)] + [sum(looped[start:start + length]) for start, length in product(range(0, L), range(1, L))]

# finds the number of common elements in given two lists, which are in this case, the created theoretical spectrum
# and the given experimental spectrum.
# this function is an essential piece of the bounding step.
def score(list1:list, list2:list):
    list3 = list(filter(lambda x: x in list1,list2))
    return len(list3)

# in a list of possible peptides, first, sorts them with respect to their scores, then returns the first N peptides.
# if there are peptides with same score, especially the Nth peptide, the function will return more than N elements.
def trim(list1:list,N:int):
    list1.sort(key=itemgetter(1), reverse=True)
    k=N-1
    while True:
        if k+1 >= len(list1):
            return list1
        elif list1[k+1][1] == list1[k][1]:
            k +=1
        else:
            break
    return list1[0:k+1]

"""
    TEH PSEUDOCODE OF THE FOLLOWING FUNCTION
LeaderboardCyclopeptideSequencing(Spectrum, N)
        Leaderboard ← set containing only the empty peptide
        LeaderPeptide ← empty peptide
        while Leaderboard is non-empty
            Leaderboard ← Expand(Leaderboard)(with every aa)
                for each Peptide in Leaderboard
                    if Mass(Peptide) = ParentMass(Spectrum)
                        if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                            LeaderPeptide ← Peptide
                    else if Mass(Peptide) < ParentMass(Spectrum)
                        add Peptide to Candidate peptides
            Leaderboard ← Trim(Candidate peptides, Spectrum, N)
        output LeaderPeptide


"""
def leaderboard_cyc_seq(N:int, experimental_spec:list):
    parent_mass = max(experimental_spec)
    # every peptide will be stored with its information(sequence, score, total mass)
    leaderboard = [('', 0, 0,)]
    leaderpeptide = ('', 0, 0)
    while leaderboard:
        candidate_peptides = []
        #branching
        for peptide in leaderboard:
            seq, score_of_pep, mass_of_pep = peptide
            for aa in amino_acids:
                candidate_peptide = seq+aa
                sub_spec_of_cand_pep = subspectrum_generator(candidate_peptide)
                mass_of_cand_pep = max(sub_spec_of_cand_pep)
                candidate_peptide_info = (candidate_peptide,score(sub_spec_of_cand_pep,experimental_spec),
                                          mass_of_cand_pep)
                if candidate_peptide_info[2] == parent_mass:
                    if candidate_peptide_info[1] > leaderpeptide[1]:
                        leaderpeptide = candidate_peptide_info
                #bounding 1
                elif candidate_peptide_info[2] < parent_mass:
                    candidate_peptides.append(candidate_peptide_info)
        #bounding 2
        leaderboard = trim(candidate_peptides, N)

    return leaderpeptide

input = open("tryme.txt","r") #this file should include the experimental spectrum of a peptide as integers
Spectrum = [int(item) for item in input.read().split(" ")]
N = 1000 #how many peptides you want to keep during trimming.
# Bigger numbers will take more time, smaller numbers might exclude a good candidate. "Choose wisely"

if __name__ == "__main__":
    start = timeit.default_timer()
    results = leaderboard_cyc_seq(N,Spectrum)
    print('-'.join(map(str, results)))
    print('')
    stop = timeit.default_timer()
    print(stop - start)


