# Coursera-Bioinformatics-Specialization
Coding challenges given at Bioinformatics Specialization on Coursera

LEADERBOARD CYCLOPEPTIDE SEQUENCING
"""""""
04.04.2020
Written by Aslı Gizem Çınar

This code finds the best sequence that can most closely create the given experimental spectrum(experimental_spec)
(which is a mass spectrometry outcome of an unknown cyclopeptide.) by using branch and bound algorithm.
The branching step adds every amino acid to a peptide.(function starts with an empty peptide.)
(accepts the aas with same mass as one. I and L, K and Q)
Bounding step trims the outcome of the branching step with respect to N.
"""""""

    THE PSEUDOCODE OF THE LEADERBOARD CYCLOPEPTIDE SEQUENCING
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

