# Rawan Abu Alkhayr

# Greedy Motifs Search 
# Greedy Motif Search Algorithm to find the most probable motifs in Dna sequences

# source: https://github.com/SinM9/Bioinformatic/blob/master/5.3%20Greedy%20Motif%20Search%20with%20pseudocounts.py
# source: https://stackoverflow.com/questions/53172461/greedy-motif-search-in-python


import random

# this func to create profile matrix from given motifs
def Create_Profile(motifs):
    k = len(motifs[0])  
    t = len(motifs)  
    profile = [[0] * k for _ in range(4)]  # profile matrix

    # count nucleotide occurrences in each column
    for motif in motifs:
        for i in range(k):
            nucleotide = motif[i]
            if nucleotide == 'A':
                profile[0][i] += 1
            elif nucleotide == 'C':
                profile[1][i] += 1
            elif nucleotide == 'G':
                profile[2][i] += 1
            elif nucleotide == 'T':
                profile[3][i] += 1

    return profile


# this func to calculate the score of motifs based on mismatches with the consensus string
def Score(motifs):
    mismatch_count = 0  
    consensus_string = ""  # store the consensus string
    profile = Create_Profile(motifs)
    
    # iterate over each position in motifs to find consensus string
    for i in range(len(profile[0])):   
        max_p = 0 
        most_frequent_nucleotide = 0 
          
        # check which nucleotide is most frequent in the current column
        for j in range(4):  
            if profile[j][i] > max_p:
                max_p = profile[j][i]
                most_frequent_nucleotide = j # store the index
        
        # append the nucleotide to the consensus string
        if most_frequent_nucleotide == 0:
            consensus_string += "A"
        elif most_frequent_nucleotide == 1:
            consensus_string += "C"
        elif most_frequent_nucleotide == 2:
            consensus_string += "G"
        elif most_frequent_nucleotide == 3:
            consensus_string += "T"
    
    # count mismatches  
    for motif in motifs:
        for i in range(len(motif)):
            if consensus_string[i] != motif[i]:
                mismatch_count += 1  
    
    return mismatch_count


# this func to find the most probable k-mer 
def profile_most_probable_kmer(text, k, profile):
    max_probability = -1  
    most_probable_kmer = text[:k]  

    # loop through all possible k-mers
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]  # current k-mer
        prob = 1  # current probability  

        # calculate the probability of the k-mer based on the profile matrix
        for j in range(k):
            nucleotide = kmer[j]  # get the nucleotide at position j
            if nucleotide == 'A':
                prob *= profile[0][j]
            elif nucleotide == 'C':
                prob *= profile[1][j]
            elif nucleotide == 'G':
                prob *= profile[2][j]
            elif nucleotide == 'T':
                prob *= profile[3][j]

        # update the maximum probability
        if prob > max_probability:
            max_probability = prob
            most_probable_kmer = kmer

    return most_probable_kmer


# greedy algorthim
def GreedyMotifSearch(Dna, k, t):
    best_motifs = []
    for i in range(t):
        best_motifs.append(Dna[i][:k])   
        
    best_score = Score(best_motifs)   

    # iterate over all k-mers in the first Dna string
    for i in range(len(Dna[0]) - k + 1):   
        motifs = [Dna[0][i:i + k]]   
        
        # build the motifs for other Dna strings
        for j in range(1, t):
            profile = Create_Profile(motifs)
            most_probable_kmer = profile_most_probable_kmer(Dna[j], k, profile)
            motifs.append(most_probable_kmer)
        
        # update the best motifs 
        current_score = Score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs   

    return best_motifs

# input for t, k, and Dna strings
t = int(input("Enter number of Dna strings: "))
k = int(input("Enter k-mer length: "))  
print("Enter the Dna strings:")
Dna = [input().strip().upper() for _ in range(t)]  


# print the result
best_motifs = GreedyMotifSearch(Dna, k, t)
print("Best Motifs:", best_motifs)
