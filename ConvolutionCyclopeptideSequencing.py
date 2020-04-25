from operator import itemgetter
from itertools import product
import timeit

amino_acids =[]
aas_and_masses ={}
masses_and_aas = {}

for i in range(57,201):
    aas_and_masses[chr(i)] = int(i)
    masses_and_aas[int(i)] = chr(i)
    amino_acids.append(chr(i))

def aas_to_mass(peptide:str):
    aa_mass_list = []
    for aa in peptide:
        aa_mass_list.append(aas_and_masses[aa])
    return aa_mass_list

def mass_to_aas(peptide_masses :list):
    aas_list = []
    for mass in peptide_masses:
        aas_list.append(masses_and_aas[mass])
    return aas_list

# creates the subspectrums of a given peptide in cyclic format.
def subspectrum_generator(peptide_str:str):
    peptide = aas_to_mass(peptide_str)
    L = len(peptide)
    looped = peptide + peptide[:L-2]
    return [0, sum(peptide)] + [sum(looped[start:start + length]) for start, length in product(range(0, L), range(1, L))]


def score(list1:list, list2:list):
    list3 = list(filter(lambda x: x in list1,list2))
    return len(list3)

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

def convolution_of_a_spectrum(spectrum:list, M:int):
    convolution = {}
    spectrum.sort()
    list2 = spectrum[0:]
    i=0
    for mass in list2:
        for x in range(0,i):
            y = mass-spectrum[x]
            if 57 <= y <= 200:
                if convolution.get(y, False):
                    convolution[y] += 1
                else:
                    convolution[y] = 1
        i += 1
    frequencies = list(convolution.items())
    frequencies.sort(key=itemgetter(1), reverse=True)
    frequencies_M_highest = trim(frequencies,M)
    M_highest_aas = [x[0] for x in frequencies_M_highest]
    return M_highest_aas

def convolution_cyc_seq(N:int,M:int, experimental_spec:list):
    parent_mass = max(experimental_spec)
    # every peptide will be stored with its information(sequence, score, total mass)
    leaderboard = [('', 0, 0,)]
    leaderpeptide = ('', 0, 0)
    frequent_aas = mass_to_aas(convolution_of_a_spectrum(experimental_spec,M))
    while leaderboard:
        candidate_peptides = []
        #branching
        for peptide in leaderboard:
            seq, score_of_pep, mass_of_pep = peptide
            for aa in frequent_aas:
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

    return aas_to_mass(leaderpeptide[0])

input = open("tryme.txt","r") #this file should include the experimental spectrum of a peptide as integers
spectrum = [int(item) for item in input.read().split(" ")]
spectrum.sort()
M= 18
N= 400

if __name__ == "__main__":
    start = timeit.default_timer()
    results = convolution_cyc_seq(N,M, spectrum)
    print('-'.join(map(str, results)))
    print('')
    stop = timeit.default_timer()
    print(stop - start)

